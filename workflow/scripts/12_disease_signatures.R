#--------------------------------
# DISEASE IMMUNE SIGNATURES
#--------------------------------

library(dplyr)
library(tidyr)

print("===== DISEASE SIGNATURES STARTED =====")

#--------------------------------
# STEP 1 — Load inputs
#--------------------------------

summary_scores <- read.csv(
  snakemake@input[["summary"]],
  row.names = 1
)

annotations <- read.csv(
  snakemake@input[["annotations"]]
)

print("Loaded summary scores")
print(dim(summary_scores))

#--------------------------------
# STEP 2 — Keep only annotated clusters
#--------------------------------

summary_scores <- summary_scores[, annotations$cluster]

# Rename clusters → biological names
colnames(summary_scores) <- annotations$program_name

#--------------------------------
# STEP 3 — Convert to long format
#--------------------------------

df_long <- summary_scores %>%
  mutate(disease = rownames(summary_scores)) %>%
  pivot_longer(
    cols = -disease,
    names_to = "program",
    values_to = "mean_score"
  )

#--------------------------------
# STEP 4 — Define states
#--------------------------------

df_long <- df_long %>%
  mutate(
    state = case_when(
      mean_score > 0.4  ~ "HIGH",
      mean_score < -0.4 ~ "LOW",
      TRUE              ~ "BASELINE"
    )
  )

#--------------------------------
# STEP 5 — Save output
#--------------------------------

out_file <- snakemake@output[["signatures"]]

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

write.csv(df_long, out_file, row.names = FALSE)

print("===== DISEASE SIGNATURES COMPLETE =====")