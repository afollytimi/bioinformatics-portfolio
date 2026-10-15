# 02 · Protein-altering germline SNVs across a clinical gene panel

## Question

Four unrelated 1000 Genomes individuals were sequenced to 30x. Across the coding exons of 20 genes commonly found on hereditary cancer and cardiac panels, how many high-confidence SNVs do they carry that alter the protein, and how do those split into missense, nonsense, stop-loss, start-loss and splice-site changes? The notebook writes both the total and the breakdown.

## Data

Alignments come from the New York Genome Center's 30x high-coverage resequencing of the 1000 Genomes samples ([ENA PRJEB31736](https://www.ebi.ac.uk/ena/browser/view/PRJEB31736); Byrska-Bishop et al., *Cell* 2022), aligned with BWA-MEM to the GRCh38 analysis set. The samples are HG00096 and HG00129 (British, GBR) and HG00268 and HG00348 (Finnish, FIN), listed in [`config/samples.tsv`](config/samples.tsv). 1000 Genomes participants gave consent for open release, and the data carry no names, phenotypes or clinical records.

The panel is defined in [`config/genes.txt`](config/genes.txt): *BRCA1, BRCA2, TP53, PALB2, ATM, CHEK2, MLH1, MSH2, MSH6, PMS2, APC, LDLR, APOB, PCSK9, MYBPC3, MYH7, KCNQ1, KCNH2, SCN5A* and *CFTR*. Regions are the Ensembl canonical coding exons of each gene, padded by 50 bp to capture splice sites.

## Method

Rather than downloading four whole-genome CRAMs (roughly 15 GB each), the pipeline streams only the panel regions straight from the archive using the CRAM index, which keeps the whole project under a few gigabytes. Coverage across the panel is logged per sample so thin regions are visible before any calls are made.

Calling is joint across all four samples with `bcftools mpileup` and `bcftools call`, requiring mapping and base quality of 20. Multiallelic sites are split and left-normalised, and each variant gets a stable `chrom_pos_ref_alt` ID so that VEP output joins back cleanly. Filtering happens in two layers. Any genotype supported by fewer than ten reads is first set to missing, then sites are kept only if they are SNVs with QUAL of at least 30 and no more than one of the four genotypes missing. Masking per genotype before applying the missingness rule is what makes the depth threshold mean something per person rather than summed across the cohort.

Ensembl VEP (release 110, offline cache) annotates the survivors, choosing one consequence per variant with its default preference for MANE Select and canonical transcripts. The notebook then takes the most severe term for each variant, counts the protein-altering classes, tabulates heterozygous and homozygous carriers per person, and reports the transition to transversion ratio as a quality check on the call set.

## Run it

```bash
conda env create -f environment.yml && conda activate germline-variants
vep_install -a cf -s homo_sapiens -y GRCh38 -c ~/.vep --CONVERT    # one-off, about 25 GB
THREADS=4 bash workflow/run_pipeline.sh                             # about 30 to 60 minutes
jupyter lab notebooks/variant_summary.ipynb
```

Outputs land in `results/summary/`: consequence counts, a per-gene table, per-sample carrier counts and a bar chart. Panel coverage is in `logs/panel_coverage.txt`.

## Results

*Pending first full run.* After the notebook has run on the real data, this section records the total number of protein-altering SNVs, the consequence breakdown, the genes contributing most variants, the Ti/Tv ratio and the mean panel coverage per sample.

## Assumptions and limits

The quality, depth and missingness thresholds are reasonable defaults, not values tuned for this data, and they are set at the top of the pipeline script so their effect is easy to test. Indels are excluded by design, which leaves out frameshifts, a major class of loss-of-function variant. Consequences are transcript-level predictions: a missense label means an amino acid changes, not that the protein is harmed, and no pathogenicity scores or ClinVar classifications are applied. Four people from two European populations cannot speak to how common any variant is. Reads were not realigned; the project relies on NYGC's published alignments.
