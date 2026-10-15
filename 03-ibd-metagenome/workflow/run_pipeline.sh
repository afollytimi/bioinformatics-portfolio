#!/usr/bin/env bash
# Phylum-level taxonomic profiling of HMP2 / IBDMDB stool metagenomes (NCBI PRJNA398089).
# Steps: sample selection, ENA download, Kraken2 classification, Bracken re-estimation at phylum level.
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"

THREADS="${THREADS:-8}"
KRAKEN_DB="${KRAKEN_DB:-reference/k2_standard_08gb}"
READ_LEN="${READ_LEN:-100}"        # must match a kmer_distrib file shipped with the database
METADATA="${METADATA:-config/hmp2_metadata.csv}"
PER_GROUP="${PER_GROUP:-12}"

mkdir -p data/raw results/kraken results/bracken logs
log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }

# 0. Preconditions
if [[ ! -f "$KRAKEN_DB/hash.k2d" ]]; then
  echo "Kraken2 database not found at $KRAKEN_DB." >&2
  echo "Download the current Standard-8 build from https://benlangmead.github.io/aws-indexes/k2 and unpack it there." >&2
  exit 1
fi
[[ -f "$KRAKEN_DB/database${READ_LEN}mers.kmer_distrib" ]] || { echo "No Bracken distribution for ${READ_LEN} bp reads in $KRAKEN_DB" >&2; exit 1; }

# 1. Sample selection (one baseline sample per participant, fixed seed)
if [[ ! -f config/samples.tsv ]]; then
  [[ -f "$METADATA" ]] || { echo "Place the HMP2 metadata CSV from https://ibdmdb.org at $METADATA" >&2; exit 1; }
  curl -fsS -o config/ena_runs.tsv \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA398089&result=read_run&fields=run_accession,sample_alias,sample_title,experiment_alias,library_strategy,fastq_ftp&format=tsv&limit=0"
  python scripts/select_samples.py "$METADATA" config/ena_runs.tsv config/samples.tsv --per-group "$PER_GROUP"
fi

# 2. Download, classify and re-estimate abundance per sample
tail -n +2 config/samples.tsv | while IFS=$'\t' read -r run participant external_id diagnosis group url1 url2; do
  if [[ -s "results/bracken/${run}.bracken" ]]; then log "${run}: done already, skipping"; continue; fi
  r1="data/raw/${run}_1.fastq.gz"; r2="data/raw/${run}_2.fastq.gz"
  [[ -s "$r1" ]] || curl -fL --retry 3 -o "$r1" "https://${url1}"
  [[ -z "$url2" || -s "$r2" ]] || curl -fL --retry 3 -o "$r2" "https://${url2}"

  log "${run}: Kraken2"
  if [[ -n "$url2" ]]; then
    kraken2 --db "$KRAKEN_DB" --threads "$THREADS" --paired "$r1" "$r2" \
      --report "results/kraken/${run}.k2report" --output - 2> "logs/${run}.kraken2.log"
  else
    kraken2 --db "$KRAKEN_DB" --threads "$THREADS" "$r1" \
      --report "results/kraken/${run}.k2report" --output - 2> "logs/${run}.kraken2.log"
  fi

  log "${run}: Bracken (phylum)"
  bracken -d "$KRAKEN_DB" -i "results/kraken/${run}.k2report" -o "results/bracken/${run}.bracken" \
    -r "$READ_LEN" -l P -t 10 > "logs/${run}.bracken.log"

  # Raw reads are large; keep only the reports unless KEEP_READS=1
  [[ "${KEEP_READS:-0}" == "1" ]] || rm -f "$r1" "$r2"
done

log "Done. Open notebooks/phylum_profiles.ipynb next."
