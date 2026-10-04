rule disease_immune_composition:
    input:
        table="results/tables/sample_immune_states.csv"
    output:
        plot="results/figures/immune_states/disease_immune_state_composition.png",
        counts="results/tables/disease_immune_state_counts.csv"
    script:
        "../scripts/19_disease_immune_composition.R"