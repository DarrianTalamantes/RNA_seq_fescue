# Objective: This code will take a gtf file and annotate it

# Extracts a fasta file from a gff or gtf file
rule bedtools:
    input:
        genome = config["genome"],
        gtf = config["scallop"]["output_file"]
    params:
        threads = config["params"]["star_mapping"]["threads"],
        genome_dir = config["directories"]["genome_idx"]
    conda:
        "../Conda_Envs/transcriptome.yaml"
    threads: 32
    output:
        expand(config["directories"]["genome_idx"] + "/" + "{file}", file=star_index_files)
    shell:

bedtools getfasta -fi genome.fasta -bed predictions.gtf -fo predicted_genes.fasta
