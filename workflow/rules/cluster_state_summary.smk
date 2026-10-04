rule cluster_state_summary:
     input:
        scores="results/tables/gsva_scores.csv",
        clusters="results/tables/immune_clusters.csv"
     output:
        summary="results/tables/cluster_state_summary.csv",
        heatmap="results/figures/immune_states/cluster_state_heatmap.png"
     script:
        "../scripts/17_cluster_state_summary.R"