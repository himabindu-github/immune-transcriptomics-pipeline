rule download_fastq:
    output:
        "data/raw/{sample}_1.fastq.gz",
        "data/raw/{sample}_2.fastq.gz"
    
    threads: config["threads"]["download"]
    shell:
        """
        fasterq-dump {wildcards.sample} --split-files -e {threads}  -O data/raw
        gzip -f data/raw/{wildcards.sample}_1.fastq
        gzip -f data/raw/{wildcards.sample}_2.fastq
        """
