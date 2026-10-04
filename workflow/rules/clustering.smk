rule gene_clustering:
    input:
        vst = "results/deseq2/vst_matrix.csv",
        genes = "results/tables/genes_var.txt",
        meta = "config/sample_metadata_gsm.tsv"
    output:
        clusters = "results/tables/gene_clusters.csv",
        heatmap = "results/figures/clustering/gene_cluster_heatmap.png"
    params:
        k = 6
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/06_gene_clustering.R"