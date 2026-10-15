"""Look up the public CRAM URL for each sample in the IGSR 30x high-coverage sequence index.

Usage: python cram_urls.py SEQUENCE_INDEX SAMPLES_TSV > sample_crams.tsv

The index is tab separated with '#'-prefixed header lines; the last header line names the columns.
Output: sample<TAB>https URL of the CRAM. Stops if any sample is missing from the index.
"""
import csv
import sys


def main() -> None:
    index_path, samples_path = sys.argv[1:3]
    header, rows = None, []
    with open(index_path) as fh:
        for line in fh:
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            rows.append(line.rstrip("\n").split("\t"))
    if not header or "SAMPLE_NAME" not in header:
        sys.exit("Could not find a SAMPLE_NAME column in the index header")

    path_col = header.index("ENA_FILE_PATH") if "ENA_FILE_PATH" in header else 0
    sample_col = header.index("SAMPLE_NAME")
    crams = {}
    for r in rows:
        if len(r) > max(path_col, sample_col) and r[path_col].endswith(".cram"):
            url = r[path_col]
            url = url.replace("ftp://", "https://") if url.startswith("ftp://") else url
            if not url.startswith("http"):
                url = "https://" + url
            crams.setdefault(r[sample_col], url)

    with open(samples_path) as fh:
        samples = [row["sample"] for row in csv.DictReader(fh, delimiter="\t")]
    missing = [s for s in samples if s not in crams]
    if missing:
        sys.exit(f"Samples not in the high-coverage index: {', '.join(missing)}")
    for s in samples:
        print(f"{s}\t{crams[s]}")


if __name__ == "__main__":
    main()
