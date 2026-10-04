library(tidyverse)

df <- read.csv(snakemake@input[[1]])

ggplot(df, aes(x = activation, y = exhaustion, color = ifn)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  labs(
    title = "CD8 Immune State Space",
    x = "Activation Score",
    y = "Exhaustion Score",
    color = "IFN"
  )

ggsave(snakemake@output[[1]], width = 6, height = 5)