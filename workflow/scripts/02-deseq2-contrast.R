library(DESeq2)
library(apeglm)

print("===== DESEQ2 COMPARISON START =====")

# Load DDS
dds_file <- snakemake@input[["dds"]]
print(paste("Loading DDS from:", dds_file))

dds <- readRDS(dds_file)
print("DDS loaded successfully")

# Params
case <- snakemake@params[["case"]]
control <- snakemake@params[["control"]]

print(paste("Running comparison:", case, "vs", control))

# Coef name
coef_name <- paste0("disease_status_", case, "_vs_", control)

print("Available coefficients:")
print(resultsNames(dds))

print(paste("Looking for coef:", coef_name))

# Run DE
if (coef_name %in% resultsNames(dds)) {

  print("Using apeglm shrinkage with coef")

  res <- lfcShrink(
    dds,
    coef = coef_name,
    type = "apeglm"
  )

} else {

  print("Coef not found → using contrast")

  res <- results(
    dds,
    contrast = c("disease_status", case, control)
  )
}

# Clean results
res <- as.data.frame(res)
res$gene <- rownames(res)

print(paste("Total genes:", nrow(res)))

res <- res[!is.na(res$padj), ]
print(paste("After removing NA padj:", nrow(res)))

res <- res[order(res$padj), ]

# Save
out_file <- snakemake@output[[1]]
print(paste("Saving results to:", out_file))

write.csv(res, out_file)

print("===== DESEQ2 COMPARISON DONE =====")
