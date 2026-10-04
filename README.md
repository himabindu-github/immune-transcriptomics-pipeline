# CD8 T Cell Transcriptomics Across Disease States

A Snakemake pipeline analyzing CD8 T cell RNA-seq across four unrelated diseases, checking along the way whether the results actually hold up.

## Aim

Reanalyze CD8-sorted RNA-seq from GEO (GSE60424), healthy controls vs. ALS, Type 1 Diabetes, sepsis, and MS (sampled before and ~24h after first interferon-beta treatment), to see whether CD8 gene expression differs by disease, and whether any consistent "immune state" pattern shows up across them.

## Data

20 CD8-sorted samples, 3-4 per group. MS_pretreatment and MS_posttreatment are the same 3 patients sampled 24h apart, not separate groups. Mostly female cohort (17 of 20).

## What was done

1. Pulled the CD8 samples out of the full GEO dataset, verified against SRA.
2. Checked data quality (`workflow/scripts/qc_contamination_check.R`) before trusting anything downstream.
3. Ran DESeq2 for each disease vs. healthy, one shared model across all 20 samples.
4. Scored samples on immune gene signatures (activation, exhaustion, cytotoxicity, interferon) two ways, a simple z-score average and GSVA.
5. Clustered the top variable genes to find data-driven gene programs, and separately clustered samples into "immune states" based on their signature scores.

The raw-read pipeline (FastQC → HISAT2 → featureCounts) was run on one sample only(not all 20), to prove it works, since GEO already provides processed counts.
## Results

**Differential expression, each disease vs. healthy:**

| Comparison | DEGs |
|---|---|
| Sepsis | 1,311 |
| MS_posttreatment | 118 |
| MS_pre vs MS_post | 101 |
| ALS | 31 |
| MS_pretreatment | 24 |
| Type1_Diabetes | 13 |

Sepsis is an outlier by 10x. That's explained by the next finding, not real disease biology.

**One sample is driving a lot of this.** Checking neutrophil markers against T-cell markers flagged 2 of 3 sepsis samples as likely contaminated with neutrophils, probably real, since very sick patients can have unusual circulating neutrophils that slip through cell sorting. One sample in particular shows up as an extreme outlier in six separate places: the marker check, the global PCA (`results/deseq2/global_PCA_plot.pdf`, PC1=47% of variance), the Sepsis heatmap and volcano plot, unsupervised gene clustering (independently enriched for "myeloid leukocyte activation," p=5×10⁻¹⁴), and the program-score PCA. It also shows up directly in `program_heatmap.png`, where Sepsis is the most extreme group on exactly the two programs tied to this contamination.

**Activation and exhaustion scores don't differ significantly by disease** (p=0.79, p=0.84 — `activation_boxplot.png`, `exhaustion_boxplot.png`).

**Two signatures ("activation," "cytotoxic") turned out to share 80% of their genes**, which is why they correlate at r=0.97, not independent biology, mostly the same genes twice.

**Sex-linked genes showed up in the ALS differential expression results**, not just in clustering where I'd already excluded them. Some of the "ALS" signal is probably a sex imbalance between the ALS and healthy sub-groups being compared.

**The 3-cluster "immune state" grouping is internally consistent but not statistically tied to disease** — cross-tabulated against 6 disease groups, most cells have 0-2 samples, too sparse for any real test. Worth noting: the 3 Sepsis samples land in three different clusters, not together, the disease itself isn't a single immune state.

## Limitations

- Small groups (3-4 per disease) limit what can be said confidently almost everywhere
- Likely neutrophil contamination in 2 of 3 sepsis samples
- Sex-linked genes filtered before clustering but not before differential expression
- Sepsis samples skew older than other groups
- Alignment pipeline only run on 1 of 20 samples (by design, see above)
- A few pipeline choices (cluster count, thresholds) were reasonable defaults, not tuned

## Next steps

- Rerun Sepsis excluding the contaminated sample
- Model the MS pre/post comparison as paired instead of independent 
- Check whether an NK-cell signal found in gene clustering is a second contamination issue
- Validate the sepsis-adjusted results against an independent dataset

## How to run it

```bash
git clone <repo-url>
cd immune-transcriptomics-pipeline
snakemake --use-conda --cores 4 -n     # dry run first
snakemake --use-conda --cores 4
```

## Repo layout

```
config/       sample metadata, contrast definitions
workflow/     Snakefile, rules, scripts, conda envs
data/         GEO count matrix
reference/    genome reference (not tracked)
resources/    curated gene signatures
results/      DE tables, figures, stats
```

## Tools

Snakemake, HISAT2, featureCounts, FastQC, DESeq2, GSVA, clusterProfiler, R (tidyverse, ggplot2)
