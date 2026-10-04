library(ggplot2)
library(dplyr)

res <- read.csv(snakemake@input[[1]])

# Remove NA
res <- res %>% filter(!is.na(padj))

# Avoid Inf
res$padj[res$padj == 0] <- 1e-300

# Add labels
res <- res %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      padj < 0.1 & log2FoldChange > 1 ~ "Up_weak",
      padj < 0.1 & log2FoldChange < -1 ~ "Down_weak",
      TRUE ~ "NS"
    )
  )

top_genes <- res %>%
  filter(padj < 0.1) %>%
  arrange(padj) %>%
  head(10)

# Plot
p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1) +
  scale_color_manual(values = c(
  "Up" = "red",
  "Down" = "blue",
  "Up_weak" = "orange",
  "Down_weak" = "skyblue",
  "NS" = "grey80"
  )) +
  theme_minimal() +
  xlim(-5, 5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.1), linetype = "dotted") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text(
    data = top_genes,
    aes(label = gene),
    size = 3,
    vjust = 1.5
  ) +
  labs(
    title = paste0(snakemake@wildcards$case, " vs ", snakemake@wildcards$control),
    x = "Log2 Fold Change",
    y = "-log10(adj p-value)"
  )

# Save
dir.create(dirname(snakemake@output[[1]]), showWarnings = FALSE, recursive = TRUE)
ggsave(snakemake@output[[1]], plot = p, width = 6, height = 5)
