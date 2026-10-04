rule sample_immune_states:
    input:
        clusters="results/tables/immune_clusters.csv",
        meta="config/sample_metadata_gsm.tsv"
    output:
        table="results/tables/sample_immune_states.csv"
    script:
        "../scripts/18_sample_immune_states.R"





