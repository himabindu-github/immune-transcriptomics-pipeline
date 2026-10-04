rule immune_cluster_analysis:
    input:
        scores="results/tables/gsva_scores.csv",
        clusters="results/tables/immune_clusters.csv",
        meta="config/sample_metadata_gsm.tsv"
    output:
        composition_plot="results/figures/immune_states/cluster_composition.png",
        boxplots="results/figures/gsva/cluster_program_boxplots.png",
        summary="results/tables/cluster_summary.csv"
    
    script:
        "../scripts/immune_cluster_analysis.R"