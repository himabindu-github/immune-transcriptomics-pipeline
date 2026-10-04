#--------------------------------
# PROGRAM SUMMARY HEATMAP
#--------------------------------

library(pheatmap)

print("===== SUMMARY HEATMAP STARTED =====")

#--------------------------------
# STEP 1 — Load summary scores
#--------------------------------

summary_scores <- read.csv(
  snakemake@input[["summary"]],
  row.names = 1
)

print("Loaded summary scores:")
print(dim(summary_scores))


#--------------------------------
# STEP 2 — Remove cluster_4 (if present)
#--------------------------------

if ("cluster_4" %in% colnames(summary_scores)) {
  summary_scores$cluster_4 <- NULL
  print("Removed cluster_4")
}


#--------------------------------
# STEP 3 — Rename clusters
#--------------------------------

colnames(summary_scores) <- c(
  "Innate_Immune_Activation",
  "Metal_Ion_Stress_Response",
  "Humoral_Immunity",
  "Cytotoxic_Immune_Response",
  "Interferon_Response"
)

print("Renamed clusters:")
print(colnames(summary_scores))


#--------------------------------
# STEP 4 — Convert to matrix
#--------------------------------

mat <- as.matrix(summary_scores)


#--------------------------------
# STEP 5 — Color palette
#--------------------------------

heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)


#--------------------------------
# STEP 6 — Save as PNG
#--------------------------------

png(
  filename = snakemake@output[["heatmap"]],
  width = 1000,
  height = 800,
  res = 120
)

pheatmap(
  mat,
  color = heat_colors,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 12,
  fontsize_col = 10,
  main = "Immune Program Activity Across Diseases"
)

dev.off()

print("===== SUMMARY HEATMAP COMPLETE =====")