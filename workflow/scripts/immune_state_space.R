library(tidyverse)

immune <- read.csv(snakemake@input[[1]], stringsAsFactors = FALSE)
gsva <- read.csv(snakemake@input[[2]], stringsAsFactors = FALSE)
print(colnames(gsva))
df <- immune %>%
  inner_join(gsva %>% select(sample, ifn), by = "sample") %>%
  rename(
    activation = activation_score,
    exhaustion = exhaustion_score
  )

write.csv(df, snakemake@output[[1]], row.names = FALSE)