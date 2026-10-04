rule cluster_immune_states:
    input:
        scores="results/tables/gsva_scores.csv",
        meta="config/sample_metadata_gsm.tsv"
    output:
        heatmap="results/figures/immune_states/immune_clusters_heatmap.png",
        clusters="results/tables/immune_clusters.csv"
    
    script:
        "../scripts/16_cluster_immune_states.R"