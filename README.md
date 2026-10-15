# Bioinformatics Portfolio

Timilehin Afolayan · Bioinformatics

Three self-contained analyses built on public sequencing data. Each one starts from raw reads in a public archive and ends in a Jupyter notebook that produces the final tables and figures, so anyone can rerun the work from accession numbers alone and check the numbers for themselves.

| Project | Question | Data | Main tools |
|---|---|---|---|
| [01 RNA-seq differential expression](01-rnaseq-airway/) | Which genes respond to dexamethasone in human airway smooth muscle cells? | GEO GSE52778 (Himes et al. 2014), 4 donors, treated vs untreated | FastQC, Trimmomatic, HISAT2, featureCounts, DESeq2 (pydeseq2) |
| [02 Germline variant analysis](02-germline-variants/) | How many high-confidence protein-altering SNVs do four unrelated individuals carry across a clinical gene panel? | 1000 Genomes 30x high-coverage, ENA PRJEB31736 | samtools, bcftools, Ensembl VEP |
| [03 Gut metagenome profiling](03-ibd-metagenome/) | Do phylum-level abundances of Bacillota (Firmicutes) and Bacteroidota differ between IBD and non-IBD participants? | HMP2 / IBDMDB, NCBI PRJNA398089 (MGnify MGYS00006120) | Kraken2, Bracken, pandas, SciPy |

## How the repository is organised

Every project folder follows the same layout. `config/` holds the sample sheet and any run settings. `workflow/run_pipeline.sh` takes reads from the archive through to a count or call table. `notebooks/` contains the analysis notebook, which reads the pipeline outputs and writes final tables and figures into `results/`. Large intermediate files (FASTQ, BAM, reference genomes, databases) are ignored by git; only the code, configuration and small summary outputs are tracked.

```
01-rnaseq-airway/
├── README.md            question, data, method, results
├── environment.yml      pinned conda environment
├── config/              sample sheet
├── workflow/            command-line pipeline
├── scripts/             small helpers called by the pipeline
├── notebooks/           analysis notebook
└── results/             outputs written by the pipeline and notebook
```

## Running a project

```bash
conda env create -f 01-rnaseq-airway/environment.yml
conda activate rnaseq-airway
bash 01-rnaseq-airway/workflow/run_pipeline.sh
jupyter lab 01-rnaseq-airway/notebooks/
```

Each project README lists its disk and memory needs. The RNA-seq and variant projects run comfortably on a workstation with 16 GB of RAM; the metagenome project needs roughly 10 GB for the 8 GB Kraken2 database.

## Testing

The notebooks are checked against small synthetic datasets with a known answer built in, so the analysis code can be validated in minutes without downloading anything:

```bash
pip install pandas numpy scipy statsmodels matplotlib pydeseq2 nbformat nbclient ipykernel
python tests/run_smoke_tests.py
```

The synthetic data exist only to exercise the code. Every reported result in the project READMEs comes from a run on the real public data.

## Status

The pipelines and notebooks are complete and pass the smoke tests. Results sections in each project README are filled in from full runs on the public data; any section still marked *pending* has not yet been run end to end.

## Licence

Code is released under the MIT licence. The datasets belong to their original authors and are used under the terms of the archives that host them; please cite the original studies listed in each project.
