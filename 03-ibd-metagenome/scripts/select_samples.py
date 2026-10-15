"""Choose one baseline stool metagenome per participant from HMP2 (IBDMDB) and write config/samples.tsv.

Usage:
  python select_samples.py HMP2_METADATA_CSV ENA_FILEREPORT_TSV OUT_TSV [--per-group 12] [--seed 2026]

HMP2_METADATA_CSV  the hmp2_metadata CSV downloaded from https://ibdmdb.org (Downloads > Metadata)
ENA_FILEREPORT_TSV run list for PRJNA398089 from the ENA portal API (the pipeline fetches it)

Selection rules, fixed before looking at any abundance data:
  * metagenomics data type only
  * earliest visit per participant, so each person contributes one sample
  * IBD = Crohn's disease (CD) or ulcerative colitis (UC); control = nonIBD
  * equal-sized random draw per group with a fixed seed so the selection is reproducible
"""
import argparse
import sys

import pandas as pd

REQUIRED = ["External ID", "Participant ID", "data_type", "diagnosis", "visit_num"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("metadata")
    ap.add_argument("filereport")
    ap.add_argument("out")
    ap.add_argument("--per-group", type=int, default=12)
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()

    meta = pd.read_csv(a.metadata, low_memory=False)
    absent = [c for c in REQUIRED if c not in meta.columns]
    if absent:
        sys.exit(f"Metadata is missing expected columns {absent}; columns present: {list(meta.columns)[:25]}")

    meta = meta[meta["data_type"] == "metagenomics"].copy()
    meta["visit_num"] = pd.to_numeric(meta["visit_num"], errors="coerce")
    meta = meta.dropna(subset=["visit_num"]).sort_values("visit_num")
    first = meta.groupby("Participant ID", as_index=False).first()

    runs = pd.read_csv(a.filereport, sep="\t", dtype=str)
    if "library_strategy" in runs.columns:
        runs = runs[runs["library_strategy"] == "WGS"]

    # HMP2 External IDs are stored either as the ENA sample alias or the sample title; use whichever matches
    best_col, best_hits = None, 0
    for col in ("sample_alias", "sample_title", "experiment_alias"):
        if col in runs.columns:
            hits = first["External ID"].isin(runs[col]).sum()
            if hits > best_hits:
                best_col, best_hits = col, hits
    if not best_col:
        sys.exit("No External ID from the metadata matched sample_alias, sample_title or experiment_alias in the ENA run list")
    print(f"Matched {best_hits} of {len(first)} participants via ENA column '{best_col}'", file=sys.stderr)

    merged = first.merge(runs, left_on="External ID", right_on=best_col, how="inner")
    merged = merged.drop_duplicates(subset="Participant ID")
    merged["group"] = merged["diagnosis"].map({"CD": "IBD", "UC": "IBD", "nonIBD": "nonIBD"})
    merged = merged.dropna(subset=["group"])

    picks = []
    for group, df in merged.groupby("group"):
        if len(df) < a.per_group:
            sys.exit(f"Only {len(df)} eligible {group} participants, fewer than --per-group {a.per_group}")
        picks.append(df.sample(n=a.per_group, random_state=a.seed))
    chosen = pd.concat(picks)

    urls = chosen["fastq_ftp"].str.split(";", expand=True)
    out = pd.DataFrame({
        "run": chosen["run_accession"].values,
        "participant": chosen["Participant ID"].values,
        "external_id": chosen["External ID"].values,
        "diagnosis": chosen["diagnosis"].values,
        "group": chosen["group"].values,
        "fastq_1": urls[0].values,
        "fastq_2": urls[1].values if urls.shape[1] > 1 else "",
    }).sort_values(["group", "run"])
    out.to_csv(a.out, sep="\t", index=False)
    print(out.groupby(["group", "diagnosis"]).size().to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
