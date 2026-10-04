#############################################
# 🧬 Immune Cluster Analysis (GSVA-based)
# -------------------------------------------
# This script:
#   1. Merges:
#        - GSVA immune scores
#        - Cluster assignments
#        - Sample metadata
#
#   2. Analyzes:
#        - Cluster composition (disease vs cluster)
#        - Immune program differences across clusters
#
#   3. Generates:
#        - Cluster composition plot
#        - Program boxplots per cluster
#
# INPUT:
#   - gsva_scores.csv
#   - immune_clusters.csv
#   - sample_metadata.tsv
#
# OUTPUT:
#   - cluster_composition.png
#   - cluster_program_boxplots.png
#   - cluster_summary.csv
#############################################

library(tidyverse)

#############################################
# 🔹 1. Inputs from Snakemake
#############################################

scores_file   <- snakemake@input[["scores"]]
clusters_file <- snakemake@input[["clusters"]]
meta_file     <- snakemake@input[["meta"]]

composition_plot_file <- snakemake@output[["composition_plot"]]
boxplot_file          <- snakemake@output[["boxplots"]]
summary_file          <- snakemake@output[["summary"]]

#############################################
# 🔹 2. Load data
#############################################

scores   <- read.csv(scores_file)
clusters <- read.csv(clusters_file)
meta     <- read.delim(meta_file)

#############################################
# 🔹 3. Merge all data
#############################################

df <- scores %>%
  left_join(clusters, by = "sample") %>%
  left_join(meta, by = c("sample" = "sample_id"))

#############################################
# 🔹 4. Basic sanity checks
#############################################

# Ensure cluster is treated as categorical
df$cluster <- as.factor(df$cluster)

#############################################
# 🔹 5. Cluster composition (Disease vs Cluster)
#############################################

# Count samples per disease per cluster
cluster_composition <- df %>%
  count(cluster, disease_status)

# Save summary table
write.csv(cluster_composition, summary_file, row.names = FALSE)

#############################################
# 🔹 6. Plot: Cluster composition
#############################################

p1 <- ggplot(df, aes(x = cluster, fill = disease_status)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  ylab("Proportion") +
  ggtitle("Disease Composition Across Immune Clusters")

ggsave(composition_plot_file, plot = p1, width = 6, height = 5)

#############################################
# 🔹 7. Program distribution per cluster
#############################################

# Convert wide → long format
df_long <- df %>%
  pivot_longer(
    cols = c(activation, cytotoxic, exhaustion, ifn, treg),
    names_to = "program",
    values_to = "score"
  )

#############################################
# 🔹 8. Plot: Program boxplots
#############################################

p2 <- ggplot(df_long, aes(x = cluster, y = score, fill = cluster)) +
  geom_boxplot() +
  facet_wrap(~ program, scales = "free_y") +
  theme_minimal() +
  ggtitle("Immune Program Activity Across Clusters")

ggsave(boxplot_file, plot = p2, width = 8, height = 6)

#############################################
# 🔹 9. Optional: Statistical testing
#############################################

# Kruskal-Wallis per program
stats_results <- df_long %>%
  group_by(program) %>%
  summarise(
    p_value = kruskal.test(score ~ cluster)$p.value
  )

# Save stats
write.csv(stats_results, sub(".csv", "_stats.csv", summary_file), row.names = FALSE)

#############################################
# ✅ DONE
#############################################