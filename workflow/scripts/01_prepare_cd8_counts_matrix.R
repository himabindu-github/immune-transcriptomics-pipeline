getwd()
setwd("/Users/himabindukumdam/Documents/Documents/immune-transcriptomics-pipeline")

getwd()
list.files("data/geo_counts")


raw_counts_file <- "data/geo_counts/GSE60424_raw_counts_GRCh38.p13_NCBI.tsv"

# Read counts file into R

raw_counts_full <- read.delim(raw_counts_file, row.names = 1)

cd8_samples_meta <- read.delim("config/sample_metadata_gsm.tsv")

colnames(raw_counts_full)
cd8_samples_meta$sample_id


# subset only cd8 samples from full matrix
cd8_counts_subset <- raw_counts_full[, cd8_samples_meta$sample_id]

# check the order
all(colnames(cd8_counts_subset) == cd8_samples_meta$sample_id)

# check missing samples
missing_samples <- setdiff(colnames(cd8_counts_subset), cd8_samples_meta$sample_id)
missing_samples

# # Verify all metadata sample IDs exist in the count matrix; stop if any are missing
if(!all(cd8_samples_meta$sample_id %in% colnames(raw_counts_full))) {
  stop("some samples in metadata are missing in sample file")}


# Save cd8_count_subset to a tab-separated file with geneids and proper formatting
write.table(cd8_counts_subset,
            file = "data/geo_counts/cd8_counts.tsv",
            sep = "\t",
            quote = F,
            row.names = TRUE,
            col.names = NA)
