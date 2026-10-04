#--------------------------------
# Load libraries
#--------------------------------

library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)

print("===== ENRICHMENT STEP STARTED =====")

#--------------------------------
# Load input
#--------------------------------

clusters_df <- read.csv(snakemake@input[["clusters"]])

print("Preview of clusters:")
print(head(clusters_df))

print("Cluster distribution:")
print(table(clusters_df$cluster))


#--------------------------------
# STEP 2 — Convert SYMBOL → ENTREZ
#--------------------------------

print("Converting gene symbols → ENTREZ IDs")

all_genes <- unique(clusters_df$gene)
print(paste("Total genes:", length(all_genes)))

gene_entrez <- mapIds(
  org.Hs.eg.db,
  keys = all_genes,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

gene_map <- data.frame(
  gene = names(gene_entrez),
  entrez = gene_entrez,
  stringsAsFactors = FALSE
)

# Remove NA mappings
gene_map <- gene_map[!is.na(gene_map$entrez), ]
print(paste("Mapped genes:", nrow(gene_map)))


#--------------------------------
# Merge mapping back to clusters
#--------------------------------

clusters_mapped <- merge(
  clusters_df,
  gene_map,
  by = "gene"
)

print("Preview after mapping:")
print(head(clusters_mapped))


#--------------------------------
# STEP 3 — GO Enrichment per cluster
#--------------------------------

print("Running GO enrichment per cluster...")

clusters <- unique(clusters_mapped$cluster)
print(paste("Total clusters:", length(clusters)))


#--------------------------------
# Loop through clusters
#--------------------------------

for (cl in clusters) {
  
  print(paste("Processing cluster:", cl))
  
  #--------------------------------
  # Get genes for this cluster
  #--------------------------------
  
  cluster_genes <- clusters_mapped %>%
    filter(cluster == cl) %>%
    pull(entrez)
  
  # Remove NA (extra safety)
  cluster_genes <- cluster_genes[!is.na(cluster_genes)]
  cluster_genes <- as.character(cluster_genes)
  
  print(paste("Genes in cluster:", length(cluster_genes)))
  
  #--------------------------------
  # Prepare output file
  #--------------------------------
  
  out_file <- snakemake@output[[paste0("cluster_", cl)]]
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  
  #--------------------------------
  # Handle empty cluster
  #--------------------------------
  
  if (length(cluster_genes) == 0) {
    
    print(paste("No valid genes for cluster", cl))
    
    write.csv(data.frame(), out_file, row.names = FALSE)
    
  } else {
    
    #--------------------------------
    # Run GO enrichment
    #--------------------------------
    
    ego <- enrichGO(
      gene = cluster_genes,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    
    #--------------------------------
    # Handle empty enrichment
    #--------------------------------
    
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
      
      print(paste("No enrichment for cluster", cl))
      
      write.csv(data.frame(), out_file, row.names = FALSE)
      
    } else {
      
      print(paste("Saving enrichment for cluster", cl))
      
      write.csv(as.data.frame(ego), out_file, row.names = FALSE)
    }
  }
  
}  # ✅ THIS WAS MISSING

#--------------------------------
# DONE
#--------------------------------

print("===== ENRICHMENT STEP COMPLETE =====")