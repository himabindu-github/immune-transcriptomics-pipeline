#############################################
# 🧬 Plot Immune States + Statistical Testing (GENERALIZED)
#############################################

library(tidyverse)

#############################################
# 🔹 1. Read inputs
#############################################

scores_file <- snakemake@input[["scores"]]
meta_file <- snakemake@input[["meta"]]

boxplot_file <- snakemake@output[["boxplot"]]
boxplot2_file <- snakemake@output[["exboxplot"]]
scatter_file <- snakemake@output[["scatter"]]
state_plot_file <- snakemake@output[["state_plot"]]

stats_file <- snakemake@output[["stats_summary"]]
pairwise_act_file <- snakemake@output[["pairwise_activation"]]
pairwise_exh_file <- snakemake@output[["pairwise_exhaustion"]]

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
# 🔹 4. Identify program columns dynamically
#############################################

program_cols <- setdiff(colnames(df), c("sample", "sample_id", "disease_status"))

#############################################
# 🔹 5. Define immune state (ONLY if core programs exist)
#############################################

if (all(c("activation", "exhaustion") %in% program_cols)) {
  
  df <- df %>%
    mutate(immune_state = case_when(
      activation > 0 & exhaustion > 0 ~ "Activated_Exhausted",
      activation > 0 & exhaustion <= 0 ~ "Activated",
      activation <= 0 & exhaustion > 0 ~ "Exhausted",
      TRUE ~ "Baseline"
    ))
}

#############################################
# 🔹 6. Create directories
#############################################

dir.create(dirname(boxplot_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(scatter_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(stats_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(state_plot_file), recursive = TRUE, showWarnings = FALSE)

#############################################
# 🔹 7. Core plots (only if activation/exhaustion exist)
#############################################

if (all(c("activation", "exhaustion") %in% program_cols)) {

  p1 <- ggplot(df, aes(x = disease_status, y = activation, fill = disease_status)) +
    geom_boxplot() +
    theme_minimal() +
    ggtitle("CD8 Activation Score")

  p2 <- ggplot(df, aes(x = disease_status, y = exhaustion, fill = disease_status)) +
    geom_boxplot() +
    theme_minimal() +
    ggtitle("CD8 Exhaustion Score")

  p3 <- ggplot(df, aes(x = activation, y = exhaustion, color = disease_status)) +
    geom_point(size = 3) +
    theme_minimal() +
    ggtitle("CD8 Immune State Map")

  ggsave(boxplot_file, plot = p1, width = 6, height = 4)
  ggsave(boxplot2_file, plot = p2, width = 6, height = 4)
  ggsave(scatter_file, plot = p3, width = 6, height = 5)

  # Immune state distribution
  if ("immune_state" %in% colnames(df)) {
    p4 <- ggplot(df, aes(x = disease_status, fill = immune_state)) +
      geom_bar(position = "fill") +
      theme_minimal() +
      ggtitle("Immune State Distribution")

    ggsave(state_plot_file, plot = p4, width = 6, height = 5)
  }
}

#############################################
# 🔹 8. GENERALIZED PROGRAM BOXPLOT (NEW 🔥)
#############################################

df_long <- df %>%
  pivot_longer(
    cols = all_of(program_cols),
    names_to = "program",
    values_to = "score"
  )

p_all <- ggplot(df_long, aes(x = disease_status, y = score, fill = disease_status)) +
  geom_boxplot() +
  facet_wrap(~program, scales = "free_y") +
  theme_minimal() +
  ggtitle("Immune Programs Across Diseases")

# Save alongside main boxplot
general_plot_file <- sub("activation_boxplot", "all_programs_boxplot", boxplot_file)
ggsave(general_plot_file, plot = p_all, width = 10, height = 6)

#############################################
# 🔹 9. Statistical Testing (ONLY core programs)
#############################################

if (all(c("activation", "exhaustion") %in% program_cols)) {

  kw_activation <- kruskal.test(activation ~ disease_status, data = df)
  kw_exhaustion <- kruskal.test(exhaustion ~ disease_status, data = df)

  pairwise_activation <- pairwise.wilcox.test(
    df$activation,
    df$disease_status,
    p.adjust.method = "BH"
  )

  pairwise_exhaustion <- pairwise.wilcox.test(
    df$exhaustion,
    df$disease_status,
    p.adjust.method = "BH"
  )

  #############################################
  # Save stats
  #############################################

  sink(stats_file)

  cat("Kruskal-Wallis Test - Activation\n")
  print(kw_activation)

  cat("\nKruskal-Wallis Test - Exhaustion\n")
  print(kw_exhaustion)

  sink()

  write.csv(as.data.frame(pairwise_activation$p.value), pairwise_act_file)
  write.csv(as.data.frame(pairwise_exhaustion$p.value), pairwise_exh_file)
}

#############################################
# ✅ DONE
#############################################