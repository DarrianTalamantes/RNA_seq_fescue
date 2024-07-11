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

rule transdecoder_gff3:
    input:
        fasta = config["bedtools"]["fasta_output"], 
        gtf = config["scallop"]["output_file"]
    conda:
        "../Conda_Envs/annotation.yaml"
    threads: 24
    output:
        gff3 = config["transdecoder"]["gff3"], # gtf_to_alignment_gff3.pl makes this
    shell:
        """
        gtf_to_alignment_gff3.pl {input.gtf} > {output.gff3}
        """
rule transdecoder_longorfs:
    input:
        gff3 = config["transdecoder"]["gff3"],
        fasta = config["bedtools"]["fasta_output"]
    conda:
        "../Conda_Envs/annotation.yaml"
    threads: 24
    params:
        output_dir = config["directories"]["transdecoder_dir"]
    output:
        long_orfs = config["transdecoder"]["long_orfs"] # longorfs makes this
    shell:
        """
        TransDecoder.LongOrfs -t {input.fasta} -O {params.output_dir}
        """

rule transdecoder_predict:
    input:
        fasta = config["bedtools"]["fasta_output"], 
        long_orfs = config["transdecoder"]["long_orfs"]
    conda:
        "../Conda_Envs/annotation.yaml"
    threads: 24
    params:
        output_dir = config["directories"]["transdecoder_dir"],
        pep_file_name = "predicted_transcripts.fasta.transdecoder.pep",
        fasta_gff3_name = "predicted_transcripts.fasta.transdecoder.gff3",
        cds_name = "predicted_transcripts.fasta.transdecoder.cds",
        bed_name = "predicted_transcripts.fasta.transdecoder.bed"
    output:
        pep_file = config["transdecoder"]["pep"], # predict makes this
        fasta_gff3 = config["transdecoder"]["fasta_gff3"] # predict makes this
    shell:
        """
        TransDecoder.Predict -t {input.fasta} -O {params.output_dir}
        mv {params.pep_file_name} {params.output_dir}
        mv {params.fasta_gff3_name} {params.output_dir}
        mv {params.cds_name} {params.output_dir}
        mv {params.bed_name} {params.output_dir}
        """

rule transdecoder_map_orfs:
    input:
        fasta_gff3 = config["transdecoder"]["fasta_gff3"],
        gff3 = config["transdecoder"]["gff3"],
        fasta = config["bedtools"]["fasta_output"]
    output:
        genome_ggf3 = config["transdecoder"]["genome_gff3"]
    threads: 24
    shell:
        """
        cdna_alignment_orf_to_genome_orf.pl \
            {input.fasta_gff3} \
            {input.gff3} \
            {input.fasta} > {output.genome_ggf3}    
        """

rule interproscon:
    input:
        pep_file = config["transdecoder"]["pep"]
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
        pep_file = config["transdecoder"]["pep"]    
    output:
        blast = config["blast"]["output"]
    params:
        db = config["blast"]["params"]["db"],
        evalue = config["blast"]["params"]["evalue"],
        outfmt = config["blast"]["params"]["outfmt"],
        num_threads = config["blast"]["params"]["num_threads"]
    threads: config["blast"]["params"]["num_threads"]
    shell:
        """
        blastp -query {input.pep} -db {params.db} -out {output.blast} -evalue {params.evalue} -outfmt {params.outfmt} -num_threads {params.num_threads}
        """
