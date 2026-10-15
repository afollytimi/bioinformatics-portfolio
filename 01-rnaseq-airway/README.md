# 01 · Dexamethasone response in airway smooth muscle (bulk RNA-seq)

## Question

In GEO series GSE52778, which genes change expression after dexamethasone treatment at a Benjamini–Hochberg adjusted p below 0.01 and an absolute log2 fold change of at least 1, once donor-to-donor variation is accounted for? The answer is a gene count plus a ranked list, both written by the notebook.

## Data

Himes and colleagues treated primary airway smooth muscle cells from four donors with 1 µM dexamethasone for 18 hours and sequenced each line with and without treatment ([PLoS One 2014, PMID 24926665](https://pubmed.ncbi.nlm.nih.gov/24926665/)). The libraries are Illumina TruSeq, unstranded, paired-end. This project uses the eight untreated and dexamethasone-only runs listed in [`config/samples.tsv`](config/samples.tsv); the albuterol arms of the series are left out. These are de-identified cell-line data released publicly through GEO and SRA, and no participant-level information is involved.

| Run | GEO sample | Donor line | Treatment |
|---|---|---|---|
| SRR1039508 | GSM1275862 | N61311 | untreated |
| SRR1039509 | GSM1275863 | N61311 | dexamethasone |
| SRR1039512 | GSM1275866 | N052611 | untreated |
| SRR1039513 | GSM1275867 | N052611 | dexamethasone |
| SRR1039516 | GSM1275870 | N080611 | untreated |
| SRR1039517 | GSM1275871 | N080611 | dexamethasone |
| SRR1039520 | GSM1275874 | N061011 | untreated |
| SRR1039521 | GSM1275875 | N061011 | dexamethasone |

## Method

Reads are pulled from ENA, checked with FastQC and trimmed with Trimmomatic (TruSeq3 adapters, sliding-window quality 15, minimum length 36). HISAT2 aligns them against GRCh38 using the prebuilt index, and featureCounts assigns read pairs to Ensembl release 110 genes in unstranded mode. A guard in the pipeline stops the run if chromosome names in the BAM and GTF disagree, since that mismatch silently produces a matrix of zeros.

Differential expression is fitted with DESeq2 through its Python port, pydeseq2, using the design `~ cell + dex`. Including the donor term makes the test paired, which is the right model here because the same four lines appear in both arms. Genes are kept when at least four samples carry ten or more fragments.

The notebook also checks itself. It looks up *CRISPLD2*, the gene the original paper highlighted, together with five well-established glucocorticoid receptor targets (*FKBP5*, *TSC22D3*, *PER1*, *ZBTB16*, *KLF15*), and warns if fewer than half come back significantly up. Failing that check would point to a labelling or annotation problem rather than biology.

## Run it

```bash
conda env create -f environment.yml && conda activate rnaseq-airway
THREADS=8 bash workflow/run_pipeline.sh          # about 2 to 3 hours, 60 GB free disk
jupyter lab notebooks/differential_expression.ipynb
```

Outputs land in `results/de/`: the full DESeq2 table, the significant subset, a one-line summary, a volcano and MA plot, and a PCA plot.

## Results

*Pending first full run.* Once the notebook has been run on the real data, this section records the number of significant genes, the split between up and down, the top ten genes by adjusted p, and how many of the six reference targets were recovered.

## Assumptions and limits

The cut-offs are conventions, so the headline count will shift if either threshold moves; the complete table is saved so other cuts can be applied without refitting. At a 1% false discovery rate, about one listed gene in a hundred is expected to be spurious. Four donors give limited power, and because all four were white men, effects that differ by sex or ancestry are invisible in this design. Cultured cells after an 18-hour exposure are a model of the airway, not the airway itself. Fold changes are unshrunk estimates, which is why mean expression is reported next to them: a large fold change on a gene with a handful of reads deserves less weight.
