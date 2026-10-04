rule plot_immune_state_2d:
    input:
        "results/tables/immune_clusters.csv"
    output:
        "results/figures/immune_states/immune_state_2d_clusters.png"
    script:
        "../scripts/plot_immune_state_2d.R"