#############################################
# 🧬 Compute Immune State Scores (FINAL)
# -------------------------------------------
# Uses Z-score normalization:
#   score = mean(z-scored expression of signature genes)
#############################################

library(tidyverse)

#############################################
# 🔹 1. Read inputs from Snakemake
#############################################

vst_file <- snakemake@input[["vst"]]
sig_file <- snakemake@input[["signatures"]]
out_file <- snakemake@output[["scores"]]

#############################################
# 🔹 2. Load data
#############################################

vst <- read.csv(vst_file, row.names = 1)
sig <- read.csv(sig_file)

#############################################
# 🔹 3. Ensure consistent gene ID format
#############################################

rownames(vst) <- as.character(rownames(vst))
sig$gene <- as.character(sig$gene)

#############################################
# 🔹 4. Z-score normalization (KEY STEP)
#############################################

# Scale each gene across samples
# Result: each gene has mean=0, sd=1
vst_z <- t(scale(t(vst)))

#############################################
# 🔹 5. Define scoring function
#############################################

compute_score <- function(genes, matrix) {
  
  genes_present <- intersect(genes, rownames(matrix))
  
  if (length(genes_present) == 0) {
    warning("No matching genes found for this signature")
    return(rep(NA, ncol(matrix)))
  }
  
  subset_matrix <- matrix[genes_present, , drop = FALSE]
  
  # Mean Z-score across genes (per sample)
  scores <- colMeans(subset_matrix)
  
  return(scores)
}

#############################################
# 🔹 6. Extract gene sets
#############################################

activation_genes <- sig %>%
  filter(signature == "activation") %>%
  pull(gene)

exhaustion_genes <- sig %>%
  filter(signature == "exhaustion") %>%
  pull(gene)

#############################################
# 🔹 7. Compute scores
#############################################

activation_score <- compute_score(activation_genes, vst_z)
exhaustion_score <- compute_score(exhaustion_genes, vst_z)

#############################################
# 🔹 8. Combine results
#############################################

scores <- data.frame(
  sample = colnames(vst),
  activation_score = activation_score,
  exhaustion_score = exhaustion_score
)

#############################################
# 🔹 9. Save output
#############################################

write.csv(scores, out_file, row.names = FALSE)

#############################################
# ✅ DONE
#############################################