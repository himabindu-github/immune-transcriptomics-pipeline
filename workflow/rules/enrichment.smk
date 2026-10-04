rule enrichment:
    input:
        clusters = "results/tables/gene_clusters.csv"
    output:
        # One file per cluster
        cluster_1 = "results/enrichment/cluster_1.csv",
        cluster_2 = "results/enrichment/cluster_2.csv",
        cluster_3 = "results/enrichment/cluster_3.csv",
        cluster_4 = "results/enrichment/cluster_4.csv",
        cluster_5 = "results/enrichment/cluster_5.csv",
        cluster_6 = "results/enrichment/cluster_6.csv"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/07_enrichment.R"