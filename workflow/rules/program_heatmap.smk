rule program_heatmap:
    input:
        summary = "results/tables/program_summary.csv"
    output:
        heatmap = "results/figures/gsva/program_heatmap.png"
    script:
        "../scripts/program_heatmap.R"