"""Execute every analysis notebook against synthetic inputs and check the planted answers.

Usage: python tests/run_smoke_tests.py
Executed copies of the notebooks are written to the temporary folder, never back into the repo.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import nbformat
import pandas as pd
from nbclient import NotebookClient

REPO = Path(__file__).resolve().parents[1]
NOTEBOOKS = {
    "01-rnaseq-airway": "differential_expression.ipynb",
    "02-germline-variants": "variant_summary.ipynb",
    "03-ibd-metagenome": "phylum_profiles.ipynb",
}


def run(project: str, project_dir: Path) -> None:
    src = REPO / project / "notebooks" / NOTEBOOKS[project]
    nb = nbformat.read(src, as_version=4)
    os.environ["PROJECT_DIR"] = str(project_dir)
    NotebookClient(nb, timeout=900, kernel_name="python3",
                   resources={"metadata": {"path": str(project_dir)}}).execute()
    nbformat.write(nb, project_dir / f"executed_{src.name}")
    print(f"  executed {src.name}")


def check(project: str, project_dir: Path) -> None:
    if project == "01-rnaseq-airway":
        sig = pd.read_csv(project_dir / "results/de/deseq2_significant.csv", index_col=0)
        planted = {f"ENSG{i:011d}" for i in range(100)}
        recall = len(planted & set(sig.index)) / len(planted)
        false_pos = len(set(sig.index) - planted)
        print(f"  recovered {recall:.0%} of planted genes, {false_pos} unplanted calls")
        assert recall >= 0.9 and false_pos <= 5
    elif project == "02-germline-variants":
        s = pd.read_csv(project_dir / "results/summary/consequence_counts.csv", index_col=0)
        total = int(s.loc["protein_altering_total", "n_variants"])
        print(f"  protein-altering SNVs counted: {total} (planted 40)")
        assert total == 40
    elif project == "03-ibd-metagenome":
        t = pd.read_csv(project_dir / "results/summary/group_comparison.csv", index_col=0)
        row = t.loc["Bacillota"]
        print(f"  Bacillota median IBD {row['median_IBD']:.3f} vs nonIBD {row['median_nonIBD']:.3f}, q = {row['q_value']:.2g}")
        assert row["median_IBD"] < row["median_nonIBD"] and row["q_value"] < 0.05


def main() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        subprocess.run([sys.executable, str(REPO / "tests/make_synthetic_data.py"), str(tmp)], check=True)
        for project in NOTEBOOKS:
            print(project)
            run(project, tmp / project)
            check(project, tmp / project)
    print("All smoke tests passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
