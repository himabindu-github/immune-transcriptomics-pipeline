#############################################
# 🧬 Correlation Between Immune Programs
#############################################

library(tidyverse)

#############################################
# 🔹 1. Inputs
#############################################

scores_file <- snakemake@input[["scores"]]
out_corr <- snakemake@output[["corr_matrix"]]
out_plot <- snakemake@output[["corr_plot"]]

#############################################
# 🔹 2. Load data
#############################################

df <- read.csv(scores_file)

#############################################
# 🔹 3. Keep only numeric program columns
#############################################

program_df <- df %>%
  select(-sample)

#############################################
# 🔹 4. Compute correlation matrix
#############################################

corr_mat <- cor(program_df, method = "spearman")

#############################################
# 🔹 5. Save matrix
#############################################

write.csv(corr_mat, out_corr)

#############################################
# 🔹 6. Plot heatmap
#############################################

corr_long <- as.data.frame(as.table(corr_mat))

p <- ggplot(corr_long, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Immune Program Correlation")

ggsave(out_plot, plot = p, width = 6, height = 5)

#############################################
# ✅ DONE
#############################################