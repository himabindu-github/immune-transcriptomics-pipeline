library(DESeq2)

counts <- read.delim(snakemake@input[["counts"]], row.names=1)
meta <- read.delim(snakemake@input[["meta"]], sep="\t", header = TRUE)

meta$sample_id <- trimws(meta$sample_id)
rownames(meta) <- meta$sample_id
meta <- meta[colnames(counts), ]

dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ disease_status
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

vsd <- vst(dds, blind=TRUE)

png(snakemake@output[[1]])
plotPCA(vsd, intgroup="disease_status")
dev.off()
