#!/usr/bin/env bash
# RNA-seq pipeline for GSE52778 (airway smooth muscle, dexamethasone vs untreated).
# Steps: ENA download, FastQC, Trimmomatic, HISAT2 alignment, featureCounts.
# Run from anywhere inside the conda environment defined in environment.yml.
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"

THREADS="${THREADS:-8}"
REF_DIR="${REF_DIR:-reference}"
ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-110}"
SAMPLES="config/samples.tsv"
HISAT2_INDEX="$REF_DIR/grch38/genome"
GTF="$REF_DIR/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gtf"
ADAPTERS="${ADAPTERS:-${CONDA_PREFIX:-}/share/trimmomatic/adapters/TruSeq3-PE-2.fa}"

mkdir -p data/raw results/qc results/trimmed results/bam results/counts logs "$REF_DIR"

log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }

# 1. Reference: prebuilt HISAT2 GRCh38 index (Ensembl chromosome names) and the matching Ensembl GTF
if [[ ! -f "${HISAT2_INDEX}.1.ht2" ]]; then
  log "Downloading prebuilt HISAT2 GRCh38 index (about 4 GB)"
  curl -fL --retry 3 -o "$REF_DIR/grch38_genome.tar.gz" https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
  tar -xzf "$REF_DIR/grch38_genome.tar.gz" -C "$REF_DIR"
  [[ -f "${HISAT2_INDEX}.1.ht2" ]] || { echo "HISAT2 index not found at ${HISAT2_INDEX} after extraction" >&2; exit 1; }
fi

if [[ ! -f "$GTF" ]]; then
  log "Downloading Ensembl release ${ENSEMBL_RELEASE} GTF"
  curl -fL --retry 3 "https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gtf.gz" \
    | gunzip -c > "$GTF"
fi

[[ -f "$ADAPTERS" ]] || { echo "Trimmomatic adapter file not found: $ADAPTERS (set ADAPTERS=...)" >&2; exit 1; }

# 2. Per-sample download, QC, trimming and alignment
tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r run gsm cell dex; do
  r1="data/raw/${run}_1.fastq.gz"
  r2="data/raw/${run}_2.fastq.gz"

  if [[ ! -s "$r1" || ! -s "$r2" ]]; then
    log "${run}: fetching FASTQ locations from ENA"
    urls=$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp&format=tsv" | tail -n +2 | cut -f2)
    IFS=';' read -r url1 url2 <<< "$urls"
    [[ -n "${url1:-}" && -n "${url2:-}" ]] || { echo "${run}: expected paired FASTQ URLs, got '${urls}'" >&2; exit 1; }
    curl -fL --retry 3 -o "$r1" "https://${url1}"
    curl -fL --retry 3 -o "$r2" "https://${url2}"
  fi

  bam="results/bam/${run}.bam"
  if [[ -s "$bam" ]]; then
    log "${run}: BAM exists, skipping"
    continue
  fi

  log "${run}: FastQC"
  fastqc -q -t 2 -o results/qc "$r1" "$r2"

  log "${run}: Trimmomatic"
  p1="results/trimmed/${run}_1.paired.fq.gz";   u1="results/trimmed/${run}_1.unpaired.fq.gz"
  p2="results/trimmed/${run}_2.paired.fq.gz";   u2="results/trimmed/${run}_2.unpaired.fq.gz"
  trimmomatic PE -threads "$THREADS" "$r1" "$r2" "$p1" "$u1" "$p2" "$u2" \
    "ILLUMINACLIP:${ADAPTERS}:2:30:10:2:True" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
    2> "logs/${run}.trimmomatic.log"

  log "${run}: HISAT2 alignment"
  hisat2 -p "$THREADS" -x "$HISAT2_INDEX" -1 "$p1" -2 "$p2" --summary-file "logs/${run}.hisat2.txt" \
    | samtools sort -@ 4 -m 1G -o "$bam" -
  samtools index "$bam"
  samtools flagstat "$bam" > "results/qc/${run}.flagstat.txt"
done

# 3. Guard against the classic chr1 versus 1 naming mismatch between BAM and GTF
first_bam=$(ls results/bam/*.bam | head -n 1)
bam_chr=$(samtools view -H "$first_bam" | awk -F'\t' '$1=="@SQ"{sub("SN:","",$2); print $2; exit}')
gtf_chr=$(grep -v '^#' "$GTF" | head -n 1 | cut -f1)
if [[ "$bam_chr" == chr* && "$gtf_chr" != chr* ]] || [[ "$bam_chr" != chr* && "$gtf_chr" == chr* ]]; then
  echo "Chromosome naming differs between BAM (${bam_chr}) and GTF (${gtf_chr}); counts would be zero." >&2
  exit 1
fi

# 4. Gene-level counts (unstranded TruSeq libraries, count read pairs as fragments)
log "featureCounts"
featureCounts -T "$THREADS" -p --countReadPairs -s 0 -t exon -g gene_id \
  -a "$GTF" -o results/counts/gene_counts.txt results/bam/*.bam \
  2> logs/featureCounts.log

# 5. Gene ID to symbol lookup for the notebook
awk -F'\t' '$3=="gene" {
    match($9, /gene_id "[^"]+"/);   id = substr($9, RSTART+9, RLENGTH-10);
    name = id;
    if (match($9, /gene_name "[^"]+"/)) name = substr($9, RSTART+11, RLENGTH-12);
    print id "\t" name
  }' "$GTF" | sort -u | { printf 'gene_id\tgene_name\n'; cat; } > results/counts/gene_names.tsv

multiqc -q -f -o results/qc results/qc logs

log "Done. Open notebooks/differential_expression.ipynb next."
