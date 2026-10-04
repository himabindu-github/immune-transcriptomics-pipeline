#--------------------------------
# Load libraries
#--------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

print("===== PROGRAM BOXPLOT STEP STARTED =====")

#--------------------------------
# STEP 1 — Load data
#--------------------------------

program_scores <- read.csv(
  snakemake@input[["scores"]],
  row.names = 1
)

meta <- read.delim(
  snakemake@input[["meta"]],
  stringsAsFactors = FALSE
)

print("Preview program_scores:")
print(head(program_scores))

print("Preview metadata:")
print(head(meta))


#--------------------------------
# STEP 2 — Fix column names
#--------------------------------
# Ensure metadata columns match your pipeline

colnames(meta)[1] <- "sample"
colnames(meta)[2] <- "condition"

#--------------------------------
# STEP 3 — Add sample column
#--------------------------------

program_scores$sample <- rownames(program_scores)

#--------------------------------
# STEP 4 — Merge with metadata
#--------------------------------

df <- merge(program_scores, meta, by = "sample")

print("Merged data preview:")
print(head(df))


#--------------------------------
# STEP 5 — Convert to long format
#--------------------------------
# From:
# sample | cluster_1 | cluster_2 ...
# To:
# sample | condition | program | score

df_long <- df %>%
  pivot_longer(
    cols = starts_with("cluster_"),
    names_to = "program",
    values_to = "score"
  )

print("Long format preview:")
print(head(df_long))


#--------------------------------
# STEP 6 — Plot boxplots
#--------------------------------

p <- ggplot(df_long, aes(x = condition, y = score, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  facet_wrap(~ program, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    legend.position = "none"
  ) +
  labs(
    title = "Program Activity Across Conditions",
    x = "Condition",
    y = "Program Score (Z-score)"
  )

print("Plot created")


#--------------------------------
# STEP 7 — Save output
#--------------------------------

out_file <- snakemake@output[["plot"]]

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

ggsave(
  out_file,
  plot = p,
  width = 12,
  height = 8
)

print("Boxplot saved")

print("===== PROGRAM BOXPLOT COMPLETE =====")