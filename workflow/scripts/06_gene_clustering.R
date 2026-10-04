print("===== GENE CLUSTERING STARTED =====")

#--------------------------------
# Load libraries
#--------------------------------
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)

#--------------------------------
# Load inputs
#--------------------------------

# VST matrix (normalized expression)
vst <- read.csv(snakemake@input[["vst"]], row.names = 1)

# Top variable genes
genes <- read.table(
  snakemake@input[["genes"]],
  stringsAsFactors = FALSE
)[,1]

# Metadata (for annotation)
meta <- read.delim(
  snakemake@input[["meta"]],
  sep = "\t",
  stringsAsFactors = FALSE
)

#--------------------------------
# Align metadata
#--------------------------------
meta$sample_id <- trimws(meta$sample_id)
rownames(meta) <- meta$sample_id

common <- intersect(colnames(vst), rownames(meta))
vst <- vst[, common]
meta <- meta[common, ]
meta <- meta[colnames(vst), ]

#--------------------------------
# STEP — Gene ID mapping
#--------------------------------


print("Mapping ENTREZ → SYMBOL")

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(vst),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Replace NA with original IDs
gene_symbols[is.na(gene_symbols)] <- rownames(vst)

# Ensure uniqueness
gene_symbols <- make.unique(gene_symbols)

# Assign
rownames(vst) <- gene_symbols

print("Mapping complete")


#--------------------------------
# Subset genes
#--------------------------------
genes <- intersect(genes, rownames(vst))
mat <- vst[genes, ]

print(paste("Genes used:", nrow(mat)))
print(paste("Samples used:", ncol(mat)))

#--------------------------------
# Remove sex genes (important)
#--------------------------------
sex_genes <- c(
  "XIST","TSIX","DDX3Y","RPS4Y1","KDM5D",
  "USP9Y","TXLNGY","UTY","PRKY","TTTY15"
)

mat <- mat[!rownames(mat) %in% sex_genes, ]

print(paste("Genes after removing sex genes:", nrow(mat)))

#--------------------------------
# Scale data (row-wise)
#--------------------------------
mat_scaled <- t(scale(t(mat)))

#--------------------------------
# Clustering
#--------------------------------

# Distance between genes
gene_dist <- dist(mat_scaled)

# Hierarchical clustering
gene_hclust <- hclust(gene_dist, method = "ward.D2")

# Number of clusters
k <- snakemake@params[["k"]]

# Assign clusters
clusters <- cutree(gene_hclust, k = k)

#--------------------------------
# Save clusters
#--------------------------------
cluster_df <- data.frame(
  gene = rownames(mat_scaled),
  cluster = clusters
)

write.csv(cluster_df, snakemake@output[["clusters"]], row.names = FALSE)

print("Clusters saved")

#--------------------------------
# Heatmap
#--------------------------------

annotation_col <- data.frame(
  condition = meta$disease_status
)
rownames(annotation_col) <- rownames(meta)

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  clustering_method = "ward.D2",
  show_rownames = FALSE,
  filename = snakemake@output[["heatmap"]]
)

print("Heatmap saved")

print("===== GENE CLUSTERING COMPLETE =====")