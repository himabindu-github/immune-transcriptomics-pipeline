#############################################
# 🧬 Cluster Immune States (GSVA-based)
#############################################

library(tidyverse)
library(pheatmap)

#############################################
# 🔹 1. Inputs
#############################################

scores_file <- snakemake@input[["scores"]]
meta_file   <- snakemake@input[["meta"]]

heatmap_file <- snakemake@output[["heatmap"]]
clusters_file <- snakemake@output[["clusters"]]

#############################################
# 🔹 2. Load data
#############################################

scores <- read.csv(scores_file)
meta <- read.delim(meta_file)

#############################################
# 🔹 3. Merge
#############################################

df <- left_join(scores, meta, by = c("sample" = "sample_id"))

#############################################
# 🔹 4. Prepare matrix for clustering
#############################################

# Remove metadata columns, keep only GSVA scores
mat <- df %>%
  select(sample, activation, cytotoxic, exhaustion, ifn, treg) %>%
  column_to_rownames("sample") %>%
  as.matrix()

#############################################
# 🔹 5. Scale (important!)
#############################################

mat_scaled <- t(scale(t(mat)))

#############################################
# 🔹 6. Hierarchical clustering
#############################################

dist_mat <- dist(mat_scaled)
hc <- hclust(dist_mat, method = "ward.D2")

#############################################
# 🔹 7. Cut into clusters
#############################################

k <- 3   # you can try 2, 3, 4 later
clusters <- cutree(hc, k = k)

cluster_df <- data.frame(
  sample = names(clusters),
  cluster = as.factor(clusters)
)

#############################################
# 🔹 8. Save clusters
#############################################

write.csv(cluster_df, clusters_file, row.names = FALSE)

#############################################
# 🔹 9. Heatmap
#############################################

annotation <- data.frame(
  Cluster = as.factor(clusters)
)

rownames(annotation) <- rownames(mat_scaled)

pheatmap(
  mat_scaled,
  annotation_row = annotation,
  show_rownames = TRUE,
  filename = heatmap_file,
  width = 6,
  height = 6
)

#############################################
# ✅ DONE
#############################################