library(tidyverse)

#############################################
# Load data
#############################################

clusters <- read.csv(snakemake@input[["clusters"]])
meta <- read.delim(snakemake@input[["meta"]])

#############################################
# Merge metadata
#############################################

df <- left_join(
  clusters,
  meta,
  by = c("sample" = "sample_id")
)

#############################################
# Add biological labels
#############################################

df <- df %>%
  mutate(
    immune_state = case_when(
      cluster == 1 ~ "Quiescent-Regulated",
      cluster == 2 ~ "Activated-Exhausted",
      cluster == 3 ~ "Effector-Inflammatory"
    )
  )

#############################################
# Select columns
#############################################

df <- df %>%
  select(
    sample,
    disease_status,
    cluster,
    immune_state
  )

#############################################
# Save
#############################################

write.csv(
  df,
  snakemake@output[["table"]],
  row.names = FALSE
)