library(tidyverse)

# -----------------------------
# Load files
# -----------------------------

state <- read.csv("results/tables/immune_state_space.csv")
clusters <- read.csv("results/tables/immune_clusters.csv")

# -----------------------------
# Merge cluster info
# -----------------------------

df <- state %>%
  left_join(clusters, by = "sample")

# Debug check
print(colnames(df))
print(head(df))

# Convert cluster to factor
df$cluster <- as.factor(df$cluster)

# -----------------------------
# Plot
# -----------------------------

p <- ggplot(df, aes(x = activation, y = exhaustion)) +

  # Quadrant lines
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +

  # Points
  geom_point(
    aes(color = ifn, shape = cluster),
    size = 4,
    alpha = 0.9
  ) +

  # IFN gradient
  scale_color_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0
  ) +

  # Quadrant labels
  annotate(
    "text",
    x = 1.2,
    y = 1.5,
    label = "Activated\nExhausted",
    fontface = "bold"
  ) +

  annotate(
    "text",
    x = 1.2,
    y = -1.5,
    label = "Effector",
    fontface = "bold"
  ) +

  annotate(
    "text",
    x = -1.2,
    y = -1.5,
    label = "Quiescent",
    fontface = "bold"
  ) +

  annotate(
    "text",
    x = -1.2,
    y = 1.5,
    label = "Dysfunctional",
    fontface = "bold"
  ) +

  labs(
    title = "CD8⁺ T Cell Immune State Space",
    subtitle = "Activation–Exhaustion landscape with IFN signaling and immune clusters",
    x = "Activation score",
    y = "Exhaustion score",
    color = "IFN",
    shape = "Cluster"
  ) +

  theme_minimal(base_size = 14)

# -----------------------------
# Save
# -----------------------------

ggsave(
  snakemake@output[[1]],
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)