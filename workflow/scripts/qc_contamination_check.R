# QC check: does the CD8-sorted fraction show contamination from other cell types?
#
# Rationale: inspecting the top variable genes (results/tables/genes_var.txt)
# showed a block of neutrophil/granulocyte genes (S100A8, S100A9, MPO, etc.)
# sitting just below the expected sex-linked genes. This script checks
# per-sample T-cell identity against neutrophil marker expression directly,
# to see whether specific samples carry contamination signal.
#
# T-cell identity markers: CD3D (915), CD3E (916), CD8A (925), CD8B (926)
# Neutrophil markers:      MPO (4353), ELANE (1991), S100A8 (6279), S100A9 (6280)

library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)

vst  <- read.csv(snakemake@input[["vst"]], row.names = 1, check.names = FALSE)
meta <- read.delim(snakemake@input[["meta"]], check.names = FALSE)

tcell_genes  <- c("915", "916", "925", "926")
neutro_genes <- c("4353", "1991", "6279", "6280")

tcell_scores  <- colMeans(vst[rownames(vst) %in% tcell_genes, ], na.rm = TRUE)
neutro_scores <- colMeans(vst[rownames(vst) %in% neutro_genes, ], na.rm = TRUE)

qc_df <- tibble(
  sample_id       = names(tcell_scores),
  tcell_score     = tcell_scores,
  neutrophil_score = neutro_scores
) %>%
  left_join(meta, by = "sample_id")

write.csv(qc_df, snakemake@output[["table"]], row.names = FALSE)

baseline <- median(qc_df$neutrophil_score)
qc_df <- qc_df %>% mutate(elevated = neutrophil_score > (baseline + 1.0))

p <- ggplot(qc_df, aes(x = tcell_score, y = neutrophil_score, color = disease_status)) +
  geom_point(size = 3) +
  geom_hline(yintercept = baseline, linetype = "dotted", color = "grey40") +
  labs(
    title = "T-cell identity vs. neutrophil marker signal per sample",
    x = "T-cell identity score (mean VST: CD3D, CD3E, CD8A, CD8B)",
    y = "Neutrophil marker score (mean VST: MPO, ELANE, S100A8, S100A9)"
  ) +
  theme_minimal()

ggsave(snakemake@output[["plot"]], p, width = 7, height = 5.5, dpi = 150)

cat("Samples flagged as elevated (neutrophil score > baseline + 1.0):\n")
print(qc_df %>% filter(elevated) %>% select(sample_id, disease_status, neutrophil_score))
