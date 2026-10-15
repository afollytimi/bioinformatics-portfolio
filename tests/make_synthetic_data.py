"""Create small synthetic inputs, laid out exactly like real pipeline outputs, with a known answer.

These files exist only to check that the notebooks run and recover a planted signal.
They are not results and are never written into the real results/ folders.

Usage: python tests/make_synthetic_data.py OUT_DIR
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RNG = np.random.default_rng(2026)
REPO = Path(__file__).resolve().parents[1]


def copy_config(project: str, dest: Path) -> None:
    (dest / "config").mkdir(parents=True, exist_ok=True)
    for f in (REPO / project / "config").glob("*"):
        (dest / "config" / f.name).write_text(f.read_text())


def rnaseq(dest: Path) -> dict:
    """Negative binomial counts for 2,000 genes; 60 planted up and 40 planted down by 4-fold."""
    copy_config("01-rnaseq-airway", dest)
    samples = pd.read_csv(dest / "config" / "samples.tsv", sep="\t")
    n_genes = 2000
    genes = [f"ENSG{i:011d}" for i in range(n_genes)]
    names = [f"GENE{i}" for i in range(n_genes)]
    known = ["CRISPLD2", "FKBP5", "TSC22D3", "PER1", "ZBTB16", "KLF15"]
    names[: len(known)] = known

    base = RNG.lognormal(mean=5, sigma=1.5, size=n_genes)
    lfc = np.zeros(n_genes)
    lfc[:60] = 2.0           # includes the six known targets
    lfc[60:100] = -2.0
    donor = {c: RNG.normal(0, 0.3, n_genes) for c in samples["cell"].unique()}

    cols = {}
    for _, s in samples.iterrows():
        mu = base * 2 ** (donor[s["cell"]] + (lfc if s["dex"] == "trt" else 0))
        disp = 0.05
        n = 1 / disp
        p = n / (n + mu)
        cols[f"results/bam/{s['run']}.bam"] = RNG.negative_binomial(n, p)

    out = dest / "results" / "counts"
    out.mkdir(parents=True, exist_ok=True)
    fc = pd.DataFrame({"Geneid": genes, "Chr": "1", "Start": 1, "End": 1000, "Strand": "+", "Length": 1000, **cols})
    with open(out / "gene_counts.txt", "w") as fh:
        fh.write("# Program:featureCounts synthetic test file\n")
        fc.to_csv(fh, sep="\t", index=False)
    pd.DataFrame({"gene_id": genes, "gene_name": names}).to_csv(out / "gene_names.tsv", sep="\t", index=False)
    return {"planted_de_genes": 100}


def variants(dest: Path) -> dict:
    """VEP tab output and genotype table with a fixed, known consequence mix."""
    copy_config("02-germline-variants", dest)
    samples = pd.read_csv(dest / "config" / "samples.tsv", sep="\t")["sample"].tolist()
    consequences = (
        ["missense_variant"] * 30
        + ["stop_gained"] * 3
        + ["splice_donor_variant"] * 2
        + ["splice_acceptor_variant&intron_variant"] * 1
        + ["missense_variant&splice_region_variant"] * 4
        + ["synonymous_variant"] * 25
        + ["intron_variant"] * 35
    )
    rows, snvs, gts = [], [], []
    ti = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    for i, csq in enumerate(consequences):
        pos = 43_044_295 + i * 37
        ref, alt = ("A", "G") if i % 3 else ("C", "A")
        vid = f"chr17_{pos}_{ref}_{alt}"
        rows.append([vid, f"chr17:{pos}", alt, csq, "BRCA1", "MODERATE", "ENST00000357654"])
        snvs.append(["chr17", pos, vid, ref, alt])
        gts.append([vid] + [RNG.choice(["0/0", "0/1", "1/1", "./."], p=[0.4, 0.4, 0.15, 0.05]) for _ in samples])

    vep_dir = dest / "results" / "vep"
    var_dir = dest / "results" / "variants"
    vep_dir.mkdir(parents=True, exist_ok=True)
    var_dir.mkdir(parents=True, exist_ok=True)
    with open(vep_dir / "annotated_snvs.tsv", "w") as fh:
        fh.write("## ENSEMBL VARIANT EFFECT PREDICTOR synthetic test file\n")
        fh.write("#Uploaded_variation\tLocation\tAllele\tConsequence\tSYMBOL\tIMPACT\tFeature\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    pd.DataFrame(snvs, columns=["CHROM", "POS", "ID", "REF", "ALT"]).to_csv(
        var_dir / "filtered_snvs.tsv", sep="\t", index=False, header=False)
    with open(var_dir / "genotypes.tsv", "w") as fh:
        fh.write("# [1]ID" + "".join(f"\t[{i + 2}]{s}:GT" for i, s in enumerate(samples)) + "\n")
        for g in gts:
            fh.write("\t".join(g) + "\n")
    n_ti = sum((r, a) in ti for _, _, _, r, a in snvs)
    return {"protein_altering": 40, "titv": n_ti / (len(snvs) - n_ti)}


def metagenome(dest: Path) -> dict:
    """Bracken phylum tables: IBD samples planted with lower Bacillota."""
    (dest / "config").mkdir(parents=True, exist_ok=True)
    runs = []
    for i in range(24):
        diagnosis = ["CD", "UC", "nonIBD", "nonIBD"][i % 4]
        runs.append({"run": f"SRR9{i:06d}", "participant": f"P{i:03d}", "diagnosis": diagnosis,
                     "group": "nonIBD" if diagnosis == "nonIBD" else "IBD",
                     "fastq_1": "", "fastq_2": ""})
    pd.DataFrame(runs).to_csv(dest / "config" / "samples.tsv", sep="\t", index=False)

    out = dest / "results" / "bracken"
    out.mkdir(parents=True, exist_ok=True)
    phyla = {"Bacillota": 1239, "Bacteroidota": 976, "Pseudomonadota": 1224,
             "Actinomycetota": 201174, "Verrucomicrobiota": 74201}
    for r in runs:
        ibd = r["group"] == "IBD"
        alpha = np.array([20 if ibd else 45, 40, 12 if ibd else 5, 4, 2], dtype=float)
        frac = RNG.dirichlet(alpha)
        reads = (frac * 2_000_000).astype(int)
        df = pd.DataFrame({
            "name": list(phyla), "taxonomy_id": list(phyla.values()), "taxonomy_lvl": "P",
            "kraken_assigned_reads": reads, "added_reads": 0, "new_est_reads": reads,
            "fraction_total_reads": reads / reads.sum(),
        })
        df.to_csv(out / f"{r['run']}.bracken", sep="\t", index=False)
    return {"expected": "Bacillota lower in IBD"}


if __name__ == "__main__":
    root = Path(sys.argv[1])
    print(rnaseq(root / "01-rnaseq-airway"))
    print(variants(root / "02-germline-variants"))
    print(metagenome(root / "03-ibd-metagenome"))
