rule program_correlation:
    input:
        scores="results/tables/gsva_scores.csv"
    output:
        corr_matrix="results/tables/program_correlation.csv",
        corr_plot="results/figures/gsva/program_correlation_heatmap.png"
    #conda:
    #    "../envs/deseq2.yaml"
    script:
        "../scripts/15_program_correlation.R"