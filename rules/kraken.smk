# Objective: This should go between Trim.smk and star.smk.

rule build_db:
    conda:
        "../Conda_Envs/kraken.yaml"
    params:
        db_dir = config["directories"]["kraken_db"],
        db_name = config["kraken"]["db_name"],
        threads = config["kraken"]["threads"]
    threads: config["kraken"]["threads"]
    output:
        db_complete = config["kraken"]["db_name"] + "/hash.k2d"
    shell:
        """
        if [ ! -d {params.db_dir} ]; then
            mkdir -p {params.db_dir}; fi

        kraken2-build --download-library PlusPFP --db {params.db_name} --threads {params.threads}
        kraken2-build --download-taxonomy --db {params.db_dir} --threads {params.threads}
        kraken2-build --build --db {params.db_dir} --threads {params.threads}
        """

# This rule preloads the db into memory so we dont have      
rule copy_db_to_memory:
    input:
        db_complete = "/dev/shm/fungi_db/hash.k2d"  
    params:
        db_mem = "/dev/shm/fungi_db",
        db_name = config["kraken"]["db_name"]
    output:
        db_in_memory = "/dev/shm/fungi_db/hash.k2d"
    shell:
    """
        # Copy the database to shared memory if it's not already there
        if [ ! -d {params.db_mem} ]; then
            cp -r {params.db_name} {params.db_mem};
        fi
    """

# This rule runs kraken
rule kraken:
    input:
        fasta_fwd = trimmed + "/{pairs}R1.fq.gz",
        fasta_rev = trimmed + "/{pairs}R2.fq.gz",
        db_complete = "/dev/shm/fungi_db/hash.k2d"  
        db_in_memory = rules.copy_db_to_memory.output.db_in_memory

    conda:
        "../Conda_Envs/kraken.yaml"
    params:
        threads = config["kraken"]["threads"],
        db_mem = "/dev/shm/fungi_db"  # Path to the copied DB in shared memory
    threads: config["kraken"]["threads"]
    output:
        krakened = config["kraken"]["classified"] + "/krakened_{pairs}.fq.gz"
    shell:
        """
        kraken2 --paired --gzip-compressed --threads {params.threads} \
                --db {params.db_mem} \
                --memory-mapping \
                --output {output.krakened} \
                {input.fasta_fwd} {input.fasta_rev}
        """

# This kraken rule uses unclassified and classified parameters to seperate stuff.
# I like this, Wallace didnt.
# rule kraken:
#     input:
#         fasta_fwd = trimmed + "/{pairs}R1.fq.gz",
#         fasta_rev = trimmed + "/{pairs}R2.fq.gz",
#         db_complete = config["kraken"]["db_name"] + "/hash.k2d"
#     conda:
#         "../Conda_Envs/kraken.yaml"
#     params:
#         cores = "5",
#         db = config["kraken"]["db_name"],
#         classified_base = config["kraken"]["classified"] + "/krakened_{pairs}",
#         unclassified_base = config["kraken"]["unclassified"] + "/krakened_{pairs}"
#     threads: 5
#     output:
#         classified1 = config["kraken"]["classified"] + "/krakened_{pairs}R1.fq.gz",
#         classified2 = config["kraken"]["classified"] + "/krakened_{pairs}R2.fq.gz",
#         unclassified1 = config["kraken"]["unclassified"] + "/krakened_{pairs}R1.fq.gz",
#         unclassified2 = config["kraken"]["unclassified"] + "/krakened_{pairs}R2.fq.gz"
#     shell:
#         """
#         kraken2 --paired --gzip-compressed --threads {params.cores} \
#                 --db {params.db} \
#                 --classified-out {params.classified_base}#.fq.gz \
#                 --unclassified-out {params.unclassified_base}#.fq.gz \
#                 {input.fasta_fwd} {input.fasta_rev}
#         """



# Example command needed
# kraken2 --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq

