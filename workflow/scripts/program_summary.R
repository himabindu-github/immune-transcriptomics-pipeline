print("===== INPUT DEBUG =====")

print("Raw snakemake@input:")
print(snakemake@input)

print("Names of inputs:")
print(names(snakemake@input))

print("Each input separately:")

for (nm in names(snakemake@input)) {
  cat("\n---", nm, "---\n")
  print(snakemake@input[[nm]])
  cat("Type:", class(snakemake@input[[nm]]), "\n")
  cat("Length:", length(snakemake@input[[nm]]), "\n")
}

print("=======================")


library(dplyr)

print("Inputs received:")
print(snakemake@input)
# Load inputs
program_scores <- read.csv(snakemake@input[[1]], row.names = 1)
metadata <- read.csv(snakemake@input[[2]], sep = "\t")



# Add sample column
program_scores$sample <- rownames(program_scores)

# Merge
merged <- program_scores %>%
  left_join(metadata, by = c("sample" = "sample_id"))

# Compute mean per disease
summary_scores <- merged %>%
  group_by(disease_status) %>%
  summarise(across(starts_with("cluster_"), mean))

# Save
write.csv(summary_scores, snakemake@output[["summary"]], row.names = FALSE)