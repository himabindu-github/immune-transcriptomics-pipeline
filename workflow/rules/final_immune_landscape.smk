rule final_immune_landscape:
    input:
        states="results/tables/sample_immune_states.csv",
        scores="results/tables/immune_state_space.csv"
    output:
        plot="results/figures/immune_states/final_immune_landscape.png"
    script:
        "../scripts/20_final_immune_landscape.R"