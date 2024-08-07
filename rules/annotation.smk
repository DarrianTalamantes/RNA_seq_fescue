# Objective: This code will take a gtf file and annotate it.
# Could not get annotations working in snakemake, output must be used with annotations_script.sh

# Extracts a fasta file from a gff or gtf file. (makes fasta file from gtf regions file)
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

rule clean_pep_file:
    input:
        pep_file = config["transdecoder"]["pep"]
    conda:
        "../Conda_Envs/annotation.yaml"
    output:
        pep_file_clean = config["transdecoder"]["pep_clean"]
    shell:
        """
        sed 's/*//g' {input.pep_file} > {output.pep_file_clean}
        """





# Transdecoder is a little difficult maybe this program could be better
# https://github.com/Gaius-Augustus/BRAKER


## This is deactivated cause, 1 it can not match positions properly and i think its useless.
# rule transdecoder_map_orfs:
#     conda:
#         "../Conda_Envs/annotation.yaml"
#     input:
#         fasta_gff3 = config["transdecoder"]["fasta_gff3"], # Should be predicted coding regions based on consenses assembly of merged bam files
#         gff3 = config["transdecoder"]["gff3"], # transcripts 
#         fasta = config["bedtools"]["fasta_output"]
#     output:
#         genome_ggf3 = config["transdecoder"]["genome_gff3"]
#     threads: 24
#     shell:
#         """
#         cdna_alignment_orf_to_genome_orf.pl \
#             {input.fasta_gff3} \
#             {input.gff3} \
#             {input.fasta} > {output.genome_ggf3}    
#        """


# This is deactivated because I can not get it to work in snakemake. You now must run annotations_script.sh

# rule eggnog_mapper:
#     input:
#         pep_file_clean = config["transdecoder"]["pep_clean"]
#     output:
#         annotations = config["eggnog_mapper"]["output"] 
#     params:
#         ann_dir = config["directories"]["annotations"],
#         num_threads = config["eggnog_mapper"]["num_threads"]
#     conda:
#         "../Conda_Envs/eggnog.yaml"
#     threads: config["eggnog_mapper"]["num_threads"]
#     shell:
#         """
#         if [ ! -d {params.ann_dir} ]; then
#         mkdir -p {params.ann_dir}; fi
#         emapper.py -i {input.pep_file_clean} --output {output.annotations} --cpu {params.num_threads} 
#         """



# Debugging
# there is no package called ‘seqLogo’  (added to annotation.yaml) fixed
# there is no package called ‘ggplot2’  (added to annotation.yaml) fixed
# transdecoder_map_orfs ran but could not map any orf to genome.  (Was not using full genome.) fixed
# Need to update python for eggnog (Sorry, Python < 3.7 is not supported)
















# Before running this rule you will need to download the interproscan databases with 
# the script interpro_set_up.sh
# Does not work with conda and cluster combined.
# rule interproscan:
#     input:
#         pep_file = config["transdecoder"]["pep"]
#     conda:
#         "../Conda_Envs/annotation.yaml"
#     threads: 24
#     params:
#         pep_file_clean = config["transdecoder"]["pep_clean"],
#         log_file = config["interproscan"]["inter_log"]
#     output:
#         interpro_tsv = config["interproscan"]["tsv_output"],
#     shell:
#         """
#         sed 's/*//g' {input.pep_file} > {params.pep_file_clean}
#         interproscan.sh -i {params.pep_file_clean} -f tsv -o {output.interpro_tsv} -T {params.log_file} -dp
#         """


# This is blastp
# rule run_blast:
#     input:
#         pep_file = config["transdecoder"]["pep"]   
#     conda:
#         "../Conda_Envs/annotation.yaml" 
#     output:
#         blast = config["blast"]["output"]
#     params:
#         db = config["blast"]["params"]["db"],
#         evalue = config["blast"]["params"]["evalue"],
#         outfmt = config["blast"]["params"]["outfmt"],
#         num_threads = config["blast"]["params"]["num_threads"],
#         blast_dir = config["directories"]["blast"]
#     threads: config["blast"]["params"]["num_threads"]
#     shell:
#         """
#         if [ ! -d {params.blast_dir} ]; then 
#             mkdir -p {params.blast_dir}; 
#         fi
#         blastp -query {input.pep_file} -db {params.db} -out {output.blast} -evalue {params.evalue} -outfmt {params.outfmt} -num_threads {params.num_threads}
#         """