#############################################
# 🧬 Compute GSVA-based Immune Scores (GENERALIZED)
#############################################

library(tidyverse)
library(GSVA)
library(BiocParallel)

#############################################
# 🔹 1. Inputs
#############################################

vst_file <- snakemake@input[["vst"]]
sig_file <- snakemake@input[["signatures"]]
out_file <- snakemake@output[["gsva_scores"]]

#############################################
# 🔹 2. Load data
#############################################

vst <- read.csv(vst_file, row.names = 1)
sig <- read.csv(sig_file)

rownames(vst) <- as.character(rownames(vst))
sig$gene <- as.character(sig$gene)

#############################################
# 🔹 3. Build gene sets
#############################################

gene_sets <- split(sig$gene, sig$signature)

#############################################
# 🔹 4. Run GSVA (FIXED API)
#############################################

params <- gsvaParam(
  expr = as.matrix(vst),
  geneSets = gene_sets,
  kcdf = "Gaussian"
)

gsva_res <- gsva(params)

#############################################
# 🔹 5. Format output (TIDY + GENERALIZED)
#############################################

# GSVA output: gene_set × sample
gsva_df <- as.data.frame(gsva_res)

# Add gene set labels
gsva_df$signature <- rownames(gsva_df)

# Convert → tidy → wide (sample × programs)
gsva_tidy <- gsva_df %>%
  pivot_longer(
    cols = -signature,
    names_to = "sample",
    values_to = "score"
  ) %>%
  pivot_wider(
    names_from = signature,
    values_from = score
  )

# Ensure sample column is first
gsva_tidy <- gsva_tidy %>%
  relocate(sample)

#############################################
# 🔹 6. Save output
#############################################

write.csv(gsva_tidy, out_file, row.names = FALSE)

#############################################
# ✅ DONE
#############################################