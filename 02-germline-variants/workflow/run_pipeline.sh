#!/usr/bin/env bash
# Germline SNV calling on a clinical gene panel from 1000 Genomes 30x high-coverage alignments (ENA PRJEB31736).
# Steps: build panel BED, stream panel reads from public CRAMs, joint calling with bcftools,
# per-genotype depth masking, site filtering, Ensembl VEP annotation.
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"

THREADS="${THREADS:-4}"
REF_DIR="${REF_DIR:-reference}"
ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-110}"
VEP_CACHE="${VEP_CACHE:-$HOME/.vep}"
IGSR="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp"
REF="$REF_DIR/GRCh38_full_analysis_set_plus_decoy_hla.fa"
SEQ_INDEX_URL="${SEQ_INDEX_URL:-$IGSR/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index}"
GTF="$REF_DIR/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gtf"

# Filters (see README for the reasoning behind each)
MIN_QUAL=30          # site quality
MIN_GT_DP=10         # genotypes below this depth are set to missing
MAX_F_MISSING=0.25   # with four samples: at most one missing genotype per site

mkdir -p "$REF_DIR" results/bam results/variants results/vep logs
log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }

# 1. Reference genome used by NYGC for the alignments (needed to decode the CRAMs)
if [[ ! -f "$REF.fai" ]]; then
  log "Downloading GRCh38 analysis set reference (about 3 GB)"
  curl -fL --retry 3 -o "$REF" "$IGSR/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
  curl -fL --retry 3 -o "$REF.fai" "$IGSR/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
fi

# 2. Panel regions: canonical coding exons of config/genes.txt, padded 50 bp for splice sites
if [[ ! -f "$GTF" ]]; then
  curl -fL --retry 3 "https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gtf.gz" \
    | gunzip -c > "$GTF"
fi
python scripts/make_regions.py config/genes.txt "$GTF" results/panel.bed --pad 50

# 3. Map sample names to public CRAM URLs
[[ -f "$REF_DIR/sequence.index" ]] || curl -fL --retry 3 -o "$REF_DIR/sequence.index" "$SEQ_INDEX_URL"
python scripts/cram_urls.py "$REF_DIR/sequence.index" config/samples.tsv > results/sample_crams.tsv

# 4. Stream only the panel regions from each remote CRAM (htslib fetches the .crai index)
while IFS=$'\t' read -r sample url; do
  out="results/bam/${sample}.panel.bam"
  if [[ -s "$out" ]]; then log "${sample}: panel BAM exists, skipping"; continue; fi
  log "${sample}: extracting panel reads"
  samtools view -@ "$THREADS" -b -M -L results/panel.bed -T "$REF" -o "$out" "$url"
  samtools index "$out"
  samtools depth -a -b results/panel.bed "$out" \
    | awk -v s="$sample" '{n++; t+=$3; if ($3>=10) c++} END {printf "%s\tmean_depth=%.1f\tpct_ge10x=%.1f\n", s, t/n, 100*c/n}' \
    >> logs/panel_coverage.txt
done < results/sample_crams.tsv

# 5. Joint calling, normalisation and stable variant IDs
log "Calling variants"
bcftools mpileup -Ou -f "$REF" -R results/panel.bed -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 --max-depth 1000 \
    results/bam/*.panel.bam \
  | bcftools call -mv -Ou \
  | bcftools norm -m -any -f "$REF" -Ou \
  | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o results/variants/raw.vcf.gz
bcftools index -t results/variants/raw.vcf.gz

# 6. Mask low-depth genotypes, then keep high-confidence SNVs
log "Filtering"
bcftools +setGT results/variants/raw.vcf.gz -Ou -- -t q -n . -i "FMT/DP<${MIN_GT_DP}" \
  | bcftools view -v snps -i "QUAL>=${MIN_QUAL} && F_MISSING<=${MAX_F_MISSING}" -Oz -o results/variants/filtered_snvs.vcf.gz
bcftools index -t results/variants/filtered_snvs.vcf.gz
bcftools stats results/variants/filtered_snvs.vcf.gz > results/variants/filtered_snvs.stats.txt

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' results/variants/filtered_snvs.vcf.gz > results/variants/filtered_snvs.tsv
bcftools query -H -f '%ID[\t%GT]\n' results/variants/filtered_snvs.vcf.gz > results/variants/genotypes.tsv

# 7. Functional annotation (one consequence per variant; VEP's default pick order prefers MANE Select, then canonical)
if [[ ! -d "$VEP_CACHE/homo_sapiens" ]]; then
  echo "VEP cache not found in $VEP_CACHE. Install it once with:" >&2
  echo "  vep_install -a cf -s homo_sapiens -y GRCh38 -c $VEP_CACHE --CONVERT   (about 25 GB)" >&2
  exit 1
fi
log "Annotating with VEP"
vep -i results/variants/filtered_snvs.vcf.gz -o results/vep/annotated_snvs.tsv --tab \
  --cache --offline --dir_cache "$VEP_CACHE" --assembly GRCh38 --fasta "$REF" \
  --pick --symbol --canonical \
  --fields "Uploaded_variation,Location,Allele,Consequence,SYMBOL,IMPACT,Feature" \
  --fork "$THREADS" --force_overwrite --no_stats

log "Done. Open notebooks/variant_summary.ipynb next."
