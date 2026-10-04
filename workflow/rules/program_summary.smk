rule program_summary:
    input:
        scores = "results/tables/program_scores.csv",
        metadata = "config/sample_metadata_gsm.tsv"
    output:
        summary = "results/tables/program_summary.csv"
    script:
        "../scripts/program_summary.R"