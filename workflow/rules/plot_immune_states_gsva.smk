rule plot_immune_states:
    input:
        scores="results/tables/gsva_scores.csv",
        meta="config/sample_metadata_gsm.tsv"
    output:
        boxplot="results/figures/gsva/activation_boxplot.png",
        exboxplot="results/figures/gsva/exhaustion_boxplot.png",
        scatter="results/figures/immune_states/immune_state_scatter.png",
        state_plot="results/figures/immune_states/immune_state_distribution.png",
        stats_summary="results/stats/immune_stats_summary.txt",
        
        pairwise_activation="results/stats/pairwise_activation.csv",
        pairwise_exhaustion="results/stats/pairwise_exhaustion.csv"
    script:
        "../scripts/14_plot_immune_states.R"