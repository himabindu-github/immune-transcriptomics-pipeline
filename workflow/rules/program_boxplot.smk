rule program_boxplot:
    input:
        scores = "results/tables/program_scores.csv",
        meta = "config/sample_metadata_gsm.tsv"
    output:
        plot = "results/figures/gsva/program_boxplots.png"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/11_program_boxplot.R"