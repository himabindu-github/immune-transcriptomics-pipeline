library(tidyverse)
library(pheatmap)

#############################################
# Load inputs
#############################################

scores <- read.csv(snakemake@input[["scores"]])
clusters <- read.csv(snakemake@input[["clusters"]])

#############################################
# Merge data
#############################################

df <- left_join(scores, clusters, by = "sample")

#############################################
# Compute cluster means
#############################################

summary_df <- df %>%
  group_by(cluster) %>%
  summarise(
    activation = mean(activation),
    cytotoxic  = mean(cytotoxic),
    exhaustion = mean(exhaustion),
    ifn        = mean(ifn),
    treg       = mean(treg)
  )

#############################################
# Save summary table
#############################################

write.csv(
  summary_df,
  snakemake@output[["summary"]],
  row.names = FALSE
)

#############################################
# Prepare heatmap matrix
#############################################

mat <- summary_df %>%
  column_to_rownames("cluster") %>%
  as.matrix()

#############################################
# Plot heatmap
#############################################

pheatmap(
  mat,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize = 12,
  border_color = "grey70",
  filename = snakemake@output[["heatmap"]],
  width = 5,
  height = 4
)