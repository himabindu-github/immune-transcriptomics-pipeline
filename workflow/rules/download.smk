rule download_fastq:
    output:
        "data/raw/{sample}.fastq"
    shell:
        """
        fasterq-dump {wildcards.sample} -O data/raw
        """
