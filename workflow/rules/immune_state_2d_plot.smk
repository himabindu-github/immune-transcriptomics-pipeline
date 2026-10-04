rule immune_state_2d_plot:
    input:
        "results/tables/immune_state_space.csv"
    output:
        "results/figures/immune_states/immune_state_2d.png"
    script:
        "../scripts/immune_state_2d.R"