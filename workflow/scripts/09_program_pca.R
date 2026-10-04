#--------------------------------
# Load libraries
#--------------------------------

library(ggplot2)
library(dplyr)

print("===== PROGRAM PCA STARTED =====")

#--------------------------------
# Load data
#--------------------------------

scores <- read.csv(snakemake@input[["scores"]], row.names = 1)
meta   <- read.csv(snakemake@input[["meta"]], sep = "\t")

print("Preview scores:")
print(head(scores))

print("Preview metadata:")
print(head(meta))


#--------------------------------
# Match sample names
#--------------------------------

# Ensure same order
scores <- scores[meta$sample_id, ]

#--------------------------------
# PCA
#--------------------------------

pca <- prcomp(scores, scale. = TRUE)

pca_df <- as.data.frame(pca$x)

# Add condition
pca_df$condition <- meta$disease_status

print("PCA complete")


#--------------------------------
# Plot
#--------------------------------

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  ggtitle("PCA of Immune Program Scores")

#--------------------------------
# Save
#--------------------------------

out_file <- snakemake@output[["ppca"]]

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

ggsave(out_file, p, width = 6, height = 5)

print("PCA plot saved")

#--------------------------------
# Extract PCA loadings
#--------------------------------

loadings <- pca$rotation

print("PCA loadings:")
print(loadings)

# Save loadings
write.csv(
  loadings,
  "results/tables/program_pca_loadings.csv"
)
print("===== PROGRAM PCA COMPLETE =====")