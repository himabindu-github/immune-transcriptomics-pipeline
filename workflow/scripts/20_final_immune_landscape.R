library(tidyverse)
library(ggrepel)

#############################################
# Load data
#############################################

states <- read.csv(snakemake@input[["states"]])
scores <- read.csv(snakemake@input[["scores"]])

#############################################
# Merge
#############################################

df <- left_join(states, scores, by = "sample")

#############################################
# Plot
#############################################

p <- ggplot(
  df,
  aes(
    x = activation,
    y = exhaustion,
    color = ifn,
    shape = immune_state
  )
) +

  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey50") +

  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey50") +

  geom_point(size = 5, alpha = 0.9) +

  geom_text_repel(
    aes(label = disease_status),
    size = 4
  ) +

  scale_color_gradient2(
    low = "navy",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +

  labs(
    title = "CD8 Immune State Landscape",
    x = "Activation Score",
    y = "Exhaustion Score",
    color = "IFN",
    shape = "Immune State"
  ) +

  theme_minimal(base_size = 14) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    )
  )

#############################################
# Save
#############################################

ggsave(
  snakemake@output[["plot"]],
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)