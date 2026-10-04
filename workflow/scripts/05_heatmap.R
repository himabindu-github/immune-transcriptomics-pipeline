# ================================
# Load libraries
# ================================
library(pheatmap)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

print("===== HEATMAP SCRIPT STARTED =====")

# ================================
# Load inputs
# ================================
vst <- read.csv(snakemake@input[["vst"]], row.names = 1)
res <- read.csv(snakemake@input[["res"]])
meta <- read.delim(snakemake@input[["meta"]], sep = "\t")

# ================================
# Clean metadata
# ================================
meta$sample_id <- trimws(meta$sample_id)
rownames(meta) <- meta$sample_id

# ================================
# Get comparison
# ================================
case <- snakemake@wildcards$case
control <- snakemake@wildcards$control

print(paste("Comparison:", case, "vs", control))

# ================================
# Subset samples
# ================================
keep <- meta$disease_status %in% c(case, control)
meta_sub <- meta[keep, ]

print(paste("Samples kept:", nrow(meta_sub)))

# ===============================
# CONVERT VST: ENTREZ → SYMBOL
# ===============================
print("Converting Entrez IDs to Gene Symbols (VST)...")

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(vst),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Replace NA with original IDs
gene_symbols[is.na(gene_symbols)] <- rownames(vst)

# 🔥 CRITICAL FIX: ensure uniqueness
gene_symbols <- make.unique(gene_symbols)

rownames(vst) <- gene_symbols
rownames(vst) <- as.character(rownames(vst))

print(paste("Genes after mapping:", nrow(vst)))

# ================================
# Subset VST matrix
# ================================
mat <- vst[, rownames(meta_sub)]

print(paste("Matrix dimensions:", nrow(mat), "genes x", ncol(mat), "samples"))

# ===============================
# PROCESS DE RESULTS
# ===============================
print("Processing DE results...")

res$gene <- as.character(res$gene)

res$symbol <- mapIds(
  org.Hs.eg.db,
  keys = res$gene,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Clean
res <- res %>%
  filter(!is.na(symbol), !is.na(padj)) %>%
  arrange(padj)

print("Top genes:")
print(head(res))

# ================================
# Select top genes
# ================================
top_genes <- res %>%
  slice_head(n = 30)

print(paste("Top genes selected:", nrow(top_genes)))

# ================================
# Match genes safely
# ================================
valid_genes <- intersect(top_genes$symbol, rownames(mat))

print(paste("Valid genes found in matrix:", length(valid_genes)))

# Subset matrix
mat <- mat[valid_genes, ]

# ================================
# Safety check
# ================================
if (nrow(mat) < 2) {
  stop("Not enough genes for heatmap")
}

# ================================
# Ensure correct sample order
# ================================
mat <- mat[, rownames(meta_sub)]

# ================================
# Scale data (row-wise)
# ================================
mat_scaled <- t(scale(t(mat)))

# ================================
# Annotation
# ================================
annotation_col <- data.frame(
  condition = meta_sub$disease_status
)

rownames(annotation_col) <- rownames(meta_sub)

# ================================
# Output path
# ================================
out_file <- snakemake@output[[1]]

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

print(paste("Saving heatmap to:", out_file))

# ================================
# Plot heatmap
# ================================
print("Generating heatmap...")

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  fontsize_row = 6,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  scale = "none",
  main = paste(case, "vs", control),
  filename = out_file
)

# ================================
# Done
# ================================
print("===== HEATMAP GENERATED SUCCESSFULLY =====")