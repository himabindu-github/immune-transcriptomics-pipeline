print("===== DEBUG INFO =====")

print("Inputs:")
print(snakemake@input)

print("Outputs:")
print(snakemake@output)

print("Params:")
print(snakemake@params)

print("Wildcards:")
print(snakemake@wildcards)


#--------------------------------
# Load libraries
#--------------------------------

library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

print("===== IMMUNE PROGRAM: GENE VARIABILITY STARTED =====")


#---------------------------------
# STEP 0 — Load inputs
#---------------------------------

# Load VST matrix
# WHAT: normalized expression values (genes × samples)
# WHY: used for downstream variability analysis
vst <- read.csv(snakemake@input[["vst"]], row.names = 1)

# Load metadata
# WHAT: sample annotations
# WHY: needed to align samples correctly
meta <- read.delim(
  snakemake@input[["meta"]],
  sep = "\t",
  stringsAsFactors = FALSE,
  header = TRUE
)

print("DEBUG: metadata preview")
print(head(meta))
print(colnames(meta))
print(dim(meta))
#----------------------------------------
# STEP 1 — Clean metadata & align samples
#----------------------------------------

# Remove extra spaces in sample IDs
# WHY: avoid mismatches like "S1 " vs "S1"
meta$sample_id <- trimws(meta$sample_id)

# Set rownames as sample IDs
# WHY: required for matching with VST columns
rownames(meta) <- meta$sample_id


# Find common samples between VST and metadata
# WHAT: intersection of sample names
# WHY: ensures both datasets refer to same samples
common_samples <- intersect(colnames(vst), rownames(meta))


# Subset VST matrix
# WHY: remove samples not present in metadata
vst <- vst[, common_samples]

# Subset metadata
# WHY: remove samples not present in VST
meta <- meta[common_samples, ]


# Reorder metadata to EXACTLY match VST column order
# WHY: critical — ensures sample-wise correctness
meta <- meta[colnames(vst), ]


# Safety check
# WHY: stop execution if mismatch exists
if (!all(colnames(vst) == rownames(meta))) {
  stop("Sample mismatch between VST and metadata!")
}

print(paste("Samples used:", ncol(vst)))


# -------------------------------
# STEP 2 — Gene ID mapping
# -------------------------------

print("Mapping gene IDs: ENTREZ → SYMBOL")

# Map ENTREZ IDs (rownames of VST) to gene symbols
# WHY: improve interpretability + required for downstream analysis
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(vst),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Replace NA mappings
# WHY: keep genes that fail mapping
gene_symbols[is.na(gene_symbols)] <- rownames(vst)

# Ensure uniqueness
# WHY: duplicate gene names break downstream steps
gene_symbols <- make.unique(gene_symbols)

# Assign mapped names
rownames(vst) <- gene_symbols

print(paste("Total genes after mapping:", nrow(vst)))


# --------------------------------
# STEP 3 — Gene variability
# --------------------------------

print("Computing gene variability...")

# Compute variance for each gene across samples
# WHAT: measure variation in expression
# WHY: high-variance genes capture biological signal
gene_var <- apply(vst, 1, var)

# Sort genes by variance (high → low)
# WHY: rank genes by importance
gene_var_sorted <- sort(gene_var, decreasing = TRUE)

# Select top variable genes
# WHAT: top 500 genes
# WHY: reduce noise, keep strong signals
top_var_genes <- names(gene_var_sorted)[1:500]

print(paste("Top variable genes selected:", length(top_var_genes)))


# --------------------------------
# STEP 4 — Save output (IMPORTANT for Snakemake)
# --------------------------------

# Save top variable genes to file
# WHAT: one gene per line
# WHY: downstream rules (clustering, enrichment) will use this
write.table(
  top_var_genes,
  snakemake@output[["var_genes"]],
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

print("Saved: top variable genes")


#--------------------------------
# DONE
#--------------------------------

print("===== GENE VARIABILITY STEP COMPLETE =====")




