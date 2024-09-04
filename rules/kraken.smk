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
            mkdir -p {params.db_dir}; 
        fi

        kraken2-build --download-taxonomy --db {params.db_dir} --threads {params.threads}
        kraken2-build --download-library PlusPFP --db {params.db_name} --threads {params.threads}
        kraken2-build --build --db {params.db_dir} --threads {params.threads}
        """
#####################
# You might have to run this to get the database
# wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20240605.tar.gz
# Then tar the file
##################

# Rule to copy the Kraken database to shared memory
rule copy_db_to_memory:
    input:
        db_complete = config["kraken"]["db_name"]
    output:
        db_in_memory = directory(config["kraken"]["db_in_memory"])
    shell:
        """
        # Copy the database to shared memory if it's not already there
        if [ ! -d {output.db_in_memory} ]; then
            cp -r {input.db_complete} {output.db_in_memory};
        fi
        """

# Rule to run Kraken using the database in memory
rule kraken:
    input:
        fasta_fwd = trimmed + "/{pairs}R1.fq.gz",
        fasta_rev = trimmed + "/{pairs}R2.fq.gz",
        db_in_memory = rules.copy_db_to_memory.output.db_in_memory
    conda:
        "../Conda_Envs/kraken.yaml"
    params:
        threads = config["kraken"]["threads"],
    threads: config["kraken"]["threads"]
    output:
        krakened = config["kraken"]["classified"] + "/krakened_{pairs}.fq.gz"
    shell:
        """
        kraken2 --paired --gzip-compressed --threads {params.threads} \
                --db {input.db_in_memory} \
                --memory-mapping \
                --output {output.krakened} \
                {input.fasta_fwd} {input.fasta_rev}
        """



#####
# must update the conda env with conda install bioconda::krakentools
####

# This rule takes the fasta files and removes the fungal reads
rule filter_reads_excluder:
    input:
        fasta_fwd = trimmed + "/{pairs}R1.fq.gz",
        fasta_rev = trimmed + "/{pairs}R2.fq.gz",
        krakened = config["kraken"]["classified"] + "/krakened_{pairs}.fq.gz"
    conda:
        "../Conda_Envs/kraken.yaml"
    params:
        taxid_fungi = "4751",
        threads = config["kraken"]["tools_threads"]
    threads: config["kraken"]["tools_threads"]
    output:
        extracted_fwd = config["kraken"]["non_fungal"] + "/{pairs}R1.fq",
        extracted_rev = config["kraken"]["non_fungal"] + "/{pairs}R2.fq"
    shell:
        """
        extract_kraken_reads.py -k {input.krakened} -s1 {input.fasta_fwd} -s2 {input.fasta_rev} \
            --exclude --include-children --taxid {params.taxid_fungi} --threads {params.threads} \
            -o {output.extracted_fwd} -o2 {output.extracted_rev}
        """

rule filter_reads_keeper:
    input:
        fasta_fwd = trimmed + "/{pairs}R1.fq.gz",
        fasta_rev = trimmed + "/{pairs}R2.fq.gz",
        krakened = config["kraken"]["classified"] + "/krakened_{pairs}.fq.gz"
    conda:
        "../Conda_Envs/kraken.yaml"
    params:
        taxid_fungi = "4751",
        threads = config["kraken"]["tools_threads"]
    threads: config["kraken"]["tools_threads"]
    output:
        extracted_fwd = config["kraken"]["fungal"] + "/{pairs}R1.fq",
        extracted_rev = config["kraken"]["fungal"] + "/{pairs}R2.fq"
    shell:
        """
        extract_kraken_reads.py -k {input.krakened} -s1 {input.fasta_fwd} -s2 {input.fasta_rev} \
            --include-children --taxid {params.taxid_fungi} --threads {params.threads} \
            -o {output.extracted_fwd} -o2 {output.extracted_rev}
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

