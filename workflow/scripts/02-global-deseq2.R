library(DESeq2)

print("===== DESEQ2 GLOBAL START =====")

# Load input
counts_file <- snakemake@input[["counts"]]
meta_file <- snakemake@input[["meta"]]

print(paste("Loading counts from:", counts_file))
print(paste("Loading metadata from:", meta_file))

counts <- read.delim(counts_file, row.names=1)
meta <- read.delim(meta_file, sep="\t")

print(paste("Counts dimensions:", dim(counts)[1], "genes x", dim(counts)[2], "samples"))
print(paste("Metadata rows:", nrow(meta)))

# Clean metadata
meta$sample_id <- trimws(meta$sample_id)
rownames(meta) <- meta$sample_id

# Align
meta <- meta[colnames(counts), ]

print("Checking alignment between counts and metadata...")

if (!all(colnames(counts) == rownames(meta))) {
  stop("Mismatch between counts and metadata")
}

print("Alignment OK")

print("Unique disease_status values:")
print(unique(meta$disease_status))

# Create DESeq2 object
print("Creating DESeq2 dataset...")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ disease_status
)

# Filter
print("Filtering low count genes...")
before <- nrow(dds)
dds <- dds[rowSums(counts(dds)) > 10, ]
after <- nrow(dds)
print(paste("Filtered genes:", before, "→", after))

# Set reference
print("Setting reference level: Healthy_Control")
dds$disease_status <- factor(dds$disease_status)
dds$disease_status <- relevel(dds$disease_status, ref="Healthy_Control")

# Run DESeq
print("Running DESeq...")
dds <- DESeq(dds)
print("DESeq completed")

# Save DDS
dds_file <- snakemake@output[["dds"]]
print(paste("Saving DDS to:", dds_file))
saveRDS(dds, dds_file)

# VST
print("Performing VST transformation...")
vsd <- vst(dds, blind=TRUE)

vst_file <- snakemake@output[["vst"]]
print(paste("Saving VST matrix to:", vst_file))
write.csv(assay(vsd), vst_file)

# PCA
pca_file <- snakemake@output[["pca"]]
print(paste("Saving PCA to:", pca_file))

pdf(pca_file)
plotPCA(vsd, intgroup="disease_status")
dev.off()

print(paste("PCA file exists:", file.exists(pca_file)))

print("===== DESEQ2 GLOBAL DONE =====")
