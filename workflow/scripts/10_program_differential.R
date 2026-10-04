#--------------------------------
# LOAD LIBRARIES
#--------------------------------

library(dplyr)
library(tidyr)

print("===== PROGRAM DIFFERENTIAL ANALYSIS STARTED =====")

#--------------------------------
# LOAD INPUT DATA
#--------------------------------

# program scores (samples × programs)
scores <- read.csv(snakemake@input[[1]], row.names = 1)

# metadata (sample_id, disease_status)
meta <- read.delim(snakemake@input[[2]])

#--------------------------------
# PREPARE DATA
#--------------------------------

# Add sample column to scores
scores$sample_id <- rownames(scores)

# Merge scores with metadata
df <- merge(scores, meta, by = "sample_id")

print("Merged data:")
print(head(df))

#--------------------------------
# CONVERT TO LONG FORMAT
#--------------------------------

# Wide → Long
# Each row = one sample + one program

df_long <- df %>%
  pivot_longer(
    cols = starts_with("cluster_"),
    names_to = "program",
    values_to = "score"
  )

print("Long format:")
print(head(df_long))

#--------------------------------
# GET CONDITIONS
#--------------------------------

conditions <- unique(df_long$disease_status)

print("Conditions found:")
print(conditions)

#--------------------------------
# GENERATE ALL PAIRWISE COMPARISONS
#--------------------------------

comparisons <- combn(conditions, 2, simplify = FALSE)

#--------------------------------
# RUN STATISTICAL TESTS
#--------------------------------

results <- data.frame()

for (comp in comparisons) {
  
  cond1 <- comp[1]
  cond2 <- comp[2]
  
  print(paste("Comparing:", cond1, "vs", cond2))
  
  for (prog in unique(df_long$program)) {
    
    # Subset for this program
    sub <- df_long %>% filter(program == prog)
    
    # Extract values for each condition
    grp1 <- sub %>% filter(disease_status == cond1) %>% pull(score)
    grp2 <- sub %>% filter(disease_status == cond2) %>% pull(score)
    
    # Skip if too few samples
    if (length(grp1) < 2 | length(grp2) < 2) {
      next
    }
    
    # Wilcoxon test (robust for small sample size)
    test <- wilcox.test(grp1, grp2)
    
    # Determine direction of change
    effect <- ifelse(mean(grp1) > mean(grp2),
                     paste0("up_in_", cond1),
                     paste0("up_in_", cond2))
    
    # Store results
    results <- rbind(results, data.frame(
      program = prog,
      condition_1 = cond1,
      condition_2 = cond2,
      p_value = test$p.value,
      effect = effect
    ))
  }
}

#--------------------------------
# MULTIPLE TEST CORRECTION
#--------------------------------

results$p_adj <- p.adjust(results$p_value, method = "BH")

print("Final results preview:")
print(head(results))

#--------------------------------
# SAVE OUTPUT
#--------------------------------

out_file <- snakemake@output[[1]]

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

write.csv(results, out_file, row.names = FALSE)

print("===== PROGRAM DIFFERENTIAL ANALYSIS COMPLETE =====")