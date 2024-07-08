# Objective: This code will take a gtf file and annotate it

# Extracts a fasta file from a gff or gtf file
rule bedtools:
    input:
        genome = config["genome"],
        gtf = config["scallop"]["output_file"]
    conda:
        "../Conda_Envs/annotation.yaml"
    threads: 24
    output:
        fasta = config["bedtools"]["fasta_output"] 
    shell:
        """
        bedtools getfasta -fi {input.genome} -bed {input.gtf} -fo {output.fasta}
        """

rule transdecoder:
    input:
        fasta = config["bedtools"]["fasta_output"], 
        gtf = config["scallop"]["output_file"]
    conda:
        "../Conda_Envs/annotation.yaml"
    threads: 24
    output:
        gff3 = config["transdecoder"]["gff3"], # gtf_to_alignment_gff3.pl makes this
        long_orfs = ["transdecoder"]["long_orfs"], # longorfs makes this
        pep_file = ["transdecoder"]["pep"], # predict makes this
        fasta_gff3 = ["transdecoder"]["fasta_gff3"] # predict makes this
    shell:
        """
        util/gtf_to_alignment_gff3.pl {input.gtf} > {output.gff3}

        TransDecoder.LongOrfs -t {input.fasta}
        TransDecoder.Predict -t {input.fasta}
        """

rule transdecoder_map_orfs:
    input:
        fasta_gff3 = ["transdecoder"]["fasta_gff3"],
        gff3 = config["transdecoder"]["gff3"],
        fasta: config["bedtools"]["fasta_output"]
    output:
        genome_ggf3 = config["transdecoder"]["genome_gff3"]
    shell:
        """
        util/cdna_alignment_orf_to_genome_orf.pl \
            {input.fasta_gff3} \
            {input.gff3} \
            {input.fasta} > {output.genome_ggf3}    
        """

rule interproscon:
    input:
        pep_file = ["transdecoder"]["pep"]
    conda:
        "../Conda_Envs/annotation.yaml"
    threads: 24
    output:
        interpro_tsv = config["interproscan"]["tsv_output"]
    shell:
        """
        interproscan.sh -i {input.pep_file} -f tsv -o {output.interpro_tsv}
        """

rule run_blast:
    input:
        pep_file = ["transdecoder"]["pep"]    
    output:
        config["blast"]["output"]
    params:
        db=config["blast"]["params"]["db"],
        evalue=config["blast"]["params"]["evalue"],
        outfmt=config["blast"]["params"]["outfmt"],
        num_threads=config["blast"]["params"]["num_threads"]
    shell:
        """
        blastp -query {input.pep} -db {params.db} -out {output} -evalue {params.evalue} -outfmt {params.outfmt} -num_threads {params.num_threads}
        """
