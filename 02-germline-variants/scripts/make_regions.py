"""Build a BED of coding exons (plus splice-site padding) for a list of genes from an Ensembl GTF.

Usage: python make_regions.py GENES_TXT ENSEMBL_GTF OUT_BED [--pad 50]

Ensembl names chromosomes 1, 2, ... while the 1000 Genomes GRCh38 reference uses chr1, chr2, ...
so a 'chr' prefix is added. Every requested gene must be found, otherwise the script stops.
"""
import argparse
import re
import sys


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("genes")
    ap.add_argument("gtf")
    ap.add_argument("out_bed")
    ap.add_argument("--pad", type=int, default=50, help="bp added to each side of every CDS block")
    a = ap.parse_args()

    wanted = {g.strip() for g in open(a.genes) if g.strip() and not g.startswith("#")}
    name_re = re.compile(r'gene_name "([^"]+)"')
    tag_re = re.compile(r'tag "Ensembl_canonical"')

    intervals, found = [], set()
    with open(a.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "CDS":
                continue
            m = name_re.search(f[8])
            if not m or m.group(1) not in wanted or not tag_re.search(f[8]):
                continue
            chrom = f[0] if f[0].startswith("chr") else f"chr{f[0]}"
            start = max(0, int(f[3]) - 1 - a.pad)  # GTF is 1-based inclusive, BED is 0-based half-open
            end = int(f[4]) + a.pad
            intervals.append((chrom, start, end, m.group(1)))
            found.add(m.group(1))

    missing = sorted(wanted - found)
    if missing:
        sys.exit(f"Genes not found as canonical CDS in the GTF: {', '.join(missing)}")

    order = {f"chr{i}": i for i in range(1, 23)} | {"chrX": 23, "chrY": 24}
    intervals.sort(key=lambda r: (order.get(r[0], 99), r[1]))

    merged = []
    for chrom, s, e, g in intervals:
        if merged and merged[-1][0] == chrom and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([chrom, s, e, g])

    with open(a.out_bed, "w") as out:
        for chrom, s, e, g in merged:
            out.write(f"{chrom}\t{s}\t{e}\t{g}\n")
    span = sum(e - s for _, s, e, _ in merged)
    print(f"{len(found)} genes, {len(merged)} intervals, {span:,} bp", file=sys.stderr)


if __name__ == "__main__":
    main()
