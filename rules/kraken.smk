# Objective: This should go between Trim.smk and star.smk.

rule build_db:
    conda:
        "../Conda_Envs/kraken.yaml"
    params:
        db_dir = config["directories"]["kraken_db"],
        db = config["kraken"]["db_name"]
    output:
        db_complete = config["kraken"]["db_name"] + "/hash.k2d"
    shell:
    """
    if [ ! -d {params.db_dir} ]; then 
        mkdir -p {params.log_dir}; fi

    kraken2-build --download-taxonomy --db {params.db}
    kraken2-build --download-library fungi --db {params.db}
    kraken2-build --build --db {params.db}

    """


rule kraken:
    input:
        fasta_fwd=trimmed + "/{pairs}R1.fq.gz",
        fasta_rev=trimmed + "/{pairs}R2.fq.gz",
        db_complete = config["kraken"]["db_name"] + "/hash.k2d"
    conda:
        "../Conda_Envs/kraken.yaml"
    params:
        cores = "5"
        db = config["kraken"]["db_name"]
        classified_base = config["kraken"]["classified"] + "/krakened_{pairs}",
        unclassified_base = config["kraken"]["unclassified"] + "/krakened_{pairs}"
    threads: 5
    output:
        classified1 = config["kraken"]["classified"] + "/krakened_{pairs}R1.fq.gz",
        classified2 = config["kraken"]["classified"] + "/krakened_{pairs}R2.fq.gz",
        unclassified1 = config["kraken"]["unclassified"] + "/krakened_{pairs}R1.fq.gz",
        unclassified2 = config["kraken"]["unclassified"] + "/krakened_{pairs}R2.fq.gz"
    shell:
        """
        kraken2 --paired -gziped-compressed --threads {params.cores} \
                --db {params.db} \
                --classified-out {params.classified_base}#.fq.gz \
                --unclassified-out {params.unclassified_base}#.fq.gz \
                {input.fasta_fwd} {input.fasta_rev}
        """

# Example command needed
# kraken2 --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq

