rule program_annotations:
    input:
        clusters="results/tables/gene_clusters.csv"
    output:
        annotations="results/tables/program_annotations.csv"
    run:
        import pandas as pd
        import sys

        print("===== VALIDATING PROGRAM ANNOTATIONS =====")

        # Load clusters
        clusters_df = pd.read_csv(input.clusters)
        clusters = sorted(clusters_df["cluster"].unique())

        print(f"Clusters found: {clusters}")

        # Load annotations
        try:
            ann = pd.read_csv(output.annotations)
        except FileNotFoundError:
            print("ERROR: program_annotations.csv not found!")
            sys.exit(1)

        # Check required columns
        required_cols = {"cluster", "program_name"}
        if not required_cols.issubset(set(ann.columns)):
            print("ERROR: Missing required columns in annotations!")
            sys.exit(1)

        # Check cluster coverage
        annotated_clusters = sorted(ann["cluster"].unique())

        print(f"Annotated clusters: {annotated_clusters}")

        missing = set(clusters) - set(annotated_clusters)

        if len(missing) > 0:
            print(f"WARNING: Missing annotations for clusters: {missing}")

        print("===== PROGRAM ANNOTATIONS VALIDATED =====")