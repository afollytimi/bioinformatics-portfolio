# 03 · Gut microbiome phyla in IBD versus non-IBD (shotgun metagenomics)

## Question

In baseline stool metagenomes from the HMP2 inflammatory bowel disease cohort, do the relative abundances of Bacillota (formerly Firmicutes) and Bacteroidota (formerly Bacteroidetes) differ between people with IBD and non-IBD controls, and does the ratio between them shift? The notebook reports group medians and means, bootstrap intervals for the difference, effect sizes and corrected p values.

## Data

The Integrative Human Microbiome Project's IBD study (IBDMDB) followed children and adults with Crohn's disease, ulcerative colitis and non-IBD controls for up to a year, collecting stool roughly every two weeks (Lloyd-Price et al., *Nature* 2019). Shotgun reads are public under [NCBI BioProject PRJNA398089](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA398089), with MGnify analyses under study [MGYS00006120](https://www.ebi.ac.uk/metagenomics/studies/MGYS00006120). Participant metadata is downloaded separately from [ibdmdb.org](https://ibdmdb.org).

Because the cohort is longitudinal, using every sample would let a few heavily sampled people dominate. [`scripts/select_samples.py`](scripts/select_samples.py) therefore takes each participant's earliest metagenome only, then draws 12 IBD and 12 non-IBD participants with a fixed random seed. The rules were set before any abundance data were examined. Participant IDs in the cohort are already de-identified codes, and no clinical details beyond diagnosis are used.

## Method

Reads are classified with Kraken2 against the Standard-8 database, and Bracken re-estimates abundance at phylum level so that reads assigned to higher ranks are redistributed rather than lost. HMP2 reads were already host-filtered by the study team. Phylum names are harmonised, since databases built before and after NCBI's 2021 renaming disagree (Firmicutes versus Bacillota, and so on).

Each phylum averaging at least 1% across samples is compared between groups with a two-sided Mann–Whitney U test, with Benjamini–Hochberg correction across the phyla tested. The Bacillota to Bacteroidota ratio is analysed on a log2 scale. Alongside p values the notebook reports a rank-biserial effect size and a 5,000-resample bootstrap interval for the difference in medians, because with 12 people per group the size and uncertainty of an effect say more than whether it crosses a threshold.

## Run it

```bash
conda env create -f environment.yml && conda activate ibd-metagenome
# one-off: unpack the current Standard-8 Kraken2 database into reference/k2_standard_08gb
#          (https://benlangmead.github.io/aws-indexes/k2)
# one-off: save the HMP2 metadata CSV from ibdmdb.org as config/hmp2_metadata.csv
THREADS=8 bash workflow/run_pipeline.sh       # roughly 2 to 4 hours for 24 samples
jupyter lab notebooks/phylum_profiles.ipynb
```

The pipeline writes `config/samples.tsv` on its first run and reuses it afterwards, so the same 24 participants are analysed every time. Raw reads are deleted after classification unless `KEEP_READS=1` is set.

## Results

*Pending first full run.* After the notebook has run on the real data, this section records the median and mean abundance of each phylum per group, the difference with its bootstrap interval, the effect size and the corrected q value, together with the log2 Bacillota to Bacteroidota ratio.

## Assumptions and limits

Relative abundances are compositional: if one phylum falls, the others rise on paper even when their absolute numbers have not changed, so this design cannot say which group actually grew or shrank. The 8 GB Kraken2 database is a capped subset of RefSeq, so unclassified reads drop out and less-studied organisms are under-represented. Bracken assumes a fixed read length of 100 bp. Age, medication, diet and disease activity all vary across HMP2 participants and none are adjusted for in a two-group comparison. A single time point per person shows association only, never cause. Pooling Crohn's disease with ulcerative colitis adds power but can blur differences between the two. Twelve people per group is enough to detect a large shift and easily too few for a modest one.
