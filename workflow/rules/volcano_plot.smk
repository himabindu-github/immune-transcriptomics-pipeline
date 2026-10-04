rule volcano_plot:
          input: "results/deseq2/{case}_vs_{control}.csv"
          output: "results/deseq2/{case}_vs_{control}_volcano.png"
          conda: "../envs/deseq2.yaml"
          script: "../scripts/04_volcano_plot.R"
