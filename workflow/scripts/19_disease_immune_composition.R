library(tidyverse)

#############################################
# Load data
#############################################

df <- read.csv(snakemake@input[["table"]])

#############################################
# Count immune states per disease
#############################################

counts <- df %>%
  count(disease_status, immune_state)

#############################################
# Save counts
#############################################

write.csv(
  counts,
  snakemake@output[["counts"]],
  row.names = FALSE
)

#############################################
# Plot
#############################################

p <- ggplot(
  counts,
  aes(
    x = disease_status,
    y = n,
    fill = immune_state
  )
) +

  geom_bar(stat = "identity") +

  labs(
    title = "Immune State Composition Across Diseases",
    x = "Disease",
    y = "Number of samples",
    fill = "Immune State"
  ) +

  theme_minimal(base_size = 14) +

  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

#############################################
# Save
#############################################

ggsave(
  snakemake@output[["plot"]],
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)