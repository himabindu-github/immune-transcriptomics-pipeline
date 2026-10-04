rule immune_state_space_table:
    input:
        immune="results/tables/immune_scores.csv",
        gsva="results/tables/gsva_scores.csv"
    output:
        "results/tables/immune_state_space.csv"
    script:
        "../scripts/immune_state_space.R"