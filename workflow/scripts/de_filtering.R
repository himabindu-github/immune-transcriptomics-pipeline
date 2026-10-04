print("===== DEBUG INFO =====")

print("Inputs:")
print(snakemake@input)

print("Outputs:")
print(snakemake@output)

print("Params:")
print(snakemake@params)

print("Wildcards:")
print(snakemake@wildcards)



print("===== DE FILTERING STEP STARTED =====")

#--------------------------------
# Load libraries
#--------------------------------

library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
#--------------------------------
# STEP 0 — Load inputs
#--------------------------------

# Load top variable genes
# WHAT: vector of gene names (SYMBOLS)
# WHY: produced from previous step
var_genes <- read.table(
  snakemake@input[["var_genes"]],
  stringsAsFactors = FALSE
)[,1]

print(paste("Variable genes loaded:", length(var_genes)))


# Load DESeq2 results
# WHAT: contains padj, log2FC, gene IDs
# WHY: used to filter statistically significant genes
res <- read.csv(snakemake@input[["res"]])

print(paste("DE results loaded:", nrow(res)))


#--------------------------------
# STEP 1 — Prepare DE genes
#--------------------------------

# Ensure gene column is character
# WHY: avoid factor issues
res$gene <- as.character(res$gene)


# Map ENTREZ → SYMBOL
# WHY: must match gene naming used in VST / previous step
res$symbol <- mapIds(
  org.Hs.eg.db,
  keys = res$gene,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)


# Remove NA mappings
# WHY: invalid gene IDs
res <- res[!is.na(res$symbol), ]


#--------------------------------
# STEP 2 — Filter significant genes
#--------------------------------

# Keep only statistically significant genes
# WHAT: padj < 0.05
# WHY: ensures high-confidence biological signal
de_genes <- res %>%
  filter(padj < 0.1) %>%
  pull(symbol)

print(paste("Significant DE genes:", length(de_genes)))


#--------------------------------
# STEP 3 — Intersect with variable genes
#--------------------------------

# Keep genes that are BOTH:
# - highly variable
# - statistically significant
# WHY: strongest biological signals
genes_var_de <- intersect(var_genes, de_genes)

print(paste("Final genes after DE filtering:", length(genes_var_de)))


#--------------------------------
# STEP 4 — Save output
#--------------------------------

# Save gene list (one gene per line)
# WHY: used in clustering and enrichment
write.table(
  de_genes,
  snakemake@output[["de_genes"]],
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  genes_var_de,
  snakemake@output[["var_de_genes"]],
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
print("Saved: DE-filtered variable genes")


#--------------------------------
# DONE
#--------------------------------

print("===== DE FILTERING STEP COMPLETE =====")