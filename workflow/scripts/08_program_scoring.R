# This script converts gene-level expression into program-level activity.
# Each "program" = a cluster of co-expressed genes (from clustering step).
#
#   Output: sample × program matrix
#    - Rows   = samples
#    - Columns = clusters (programs)
#    - Values = activity score (Z-score based)
#
# WHY THIS IS IMPORTANT:
# Instead of analyzing thousands of genes, we summarize into a few
# biologically meaningful programs (e.g. neutrophils, cytotoxicity).

#--------------------------------
# Load libraries
#--------------------------------

library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

print("===== PROGRAM SCORING STARTED =====")

#--------------------------------
# STEP 0 — Load inputs
#--------------------------------
# INPUT 1: VST matrix
# - File: results/deseq2/vst_matrix.csv
# - Format: genes × samples
# - Values: normalized expression (variance-stabilized)
#
# Example:
#        Sample1 Sample2 Sample3
# GeneA   5.2     6.1     4.9
# GeneB   3.1     2.8     3.5
# VST matrix (genes × samples)
vst <- read.csv(snakemake@input[["vst"]], row.names = 1)
print(head(vst))

#  INPUT 2: Gene clusters
# - File: results/tables/gene_clusters.csv
# - Format: gene → cluster assignment
#
# Example:
# gene     cluster
# S100A8   1
# GZMB     5
clusters_df <- read.csv(snakemake@input[["clusters"]])

print(paste("Genes in VST:", nrow(vst)))
print(paste("Samples:", ncol(vst)))

print("Preview clusters:")
print(head(clusters_df))


#--------------------------------
# STEP 0.5 — Map ENTREZ → SYMBOL
#--------------------------------

#  WHAT:
# Convert gene identifiers from ENTREZ IDs (numeric) to gene SYMBOLS
#
# WHY:
# - VST matrix uses ENTREZ IDs
# - clustering results use gene SYMBOLS
# - mismatch prevents overlap and breaks downstream analysis
#
# OUTPUT:
# - vst matrix with gene SYMBOLS as rownames

print("Mapping ENTREZ IDs to gene symbols...")

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(vst),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Remove genes without symbol
valid_idx <- !is.na(gene_symbols)
vst <- vst[valid_idx, ]
gene_symbols <- gene_symbols[valid_idx]

# Remove duplicate mappings (keep first occurrence)
dup_idx <- duplicated(gene_symbols)
vst <- vst[!dup_idx, ]
gene_symbols <- gene_symbols[!dup_idx]

# Assign SYMBOL as rownames
rownames(vst) <- gene_symbols

print(paste("Genes after mapping:", nrow(vst)))


#--------------------------------
# STEP 1 — Match genes
#--------------------------------

#  WHAT:
# Ensure both VST matrix and cluster file use the same genes

# WHY:
# - clustering step used subset of genes
# - VST contains all genes
# - mismatch will cause errors or incorrect scoring

#  OUTPUT:
# - filtered VST matrix
# - filtered cluster dataframe

# Keep only genes present in BOTH
common_genes <- intersect(rownames(vst), clusters_df$gene)

print(paste("Common genes:", length(common_genes)))

# Subset vst to only clustered genes
vst <- vst[common_genes, ]

# # Subset clusters to valid genes
clusters_df <- clusters_df %>%
  filter(gene %in% common_genes)


#--------------------------------
# STEP 2 — Z-score normalization (per gene)
#--------------------------------

print("Scaling genes (Z-score)...")

#  WHAT:
# Standardize each gene across samples:
#   mean = 0
#   sd   = 1

# WHY:
# - genes have different expression scales
# - prevents highly expressed genes from dominating
# - focuses on relative patterns across samples

# OUTPUT:
# vst_scaled (same shape as vst, but normalized)


# WHAT: normalize each gene across samples
# WHY: remove magnitude differences, keep patterns
vst_scaled <- t(scale(t(vst)))

# Handle constant genes (sd = 0 → NA)
# Replace NA with 0 (neutral contribution)
vst_scaled[is.na(vst_scaled)] <- 0


#--------------------------------
# STEP 3 — Compute program scores
#--------------------------------

print("Computing program scores...")

#  WHAT:
# For each cluster:
#   - take all genes in that cluster
#   - compute mean expression per sample

# WHY:
# - summarizes gene set into single score
# - reduces dimensionality (genes → programs)
# - enables system-level analysis

#  OUTPUT:
# program_scores (samples × clusters)

clusters <- sort(unique(clusters_df$cluster))

print(paste("Total clusters:", length(clusters)))

# Initialize empty dataframe
# Rows = samples
program_scores <- data.frame(row.names = colnames(vst_scaled))


# Loop through clusters
for (cl in clusters) {
  
  print(paste("Processing cluster:", cl))
  
  # Get genes for this cluster
  genes_cl <- clusters_df %>%
    filter(cluster == cl) %>%
    pull(gene)
  
  # Subset matrix
  mat_cl <- vst_scaled[genes_cl, , drop = FALSE]
  
  # Safety check
  if (nrow(mat_cl) < 5) {
    print(paste("Skipping cluster", cl, "- too few genes"))
    next
  }
  
  #  WHAT:
  # mean Z-score per sample
  
  # OUTPUT:
  # vector of length = number of samples
  scores <- colMeans(mat_cl)
  
  # Add as column
  program_scores[[paste0("cluster_", cl)]] <- scores
}

print("Program scoring complete")


#--------------------------------
# STEP 4 — Save output
#--------------------------------

# OUTPUT FILE:
# results/tables/program_scores.csv

# Format:
# rows   = samples
# columns = clusters
#
# Example:
# sample     cluster_1 cluster_2 cluster_3
# GSM1       1.2       -0.3      0.8

out_file <- snakemake@output[["scores"]]

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

write.csv(
  program_scores,
  out_file,
  row.names = TRUE
)

print("Saved program scores")


#--------------------------------
# DONE
#--------------------------------

print("===== PROGRAM SCORING COMPLETE =====")