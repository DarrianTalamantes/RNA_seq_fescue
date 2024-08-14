# setting variables
genome="/scratch/drt83172/Wallace_lab/RNA_SEQ/Genome/tall_fescuev0.1.fa"
small_gtf_dir="/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/small_gtfs"
bedtools_dir="/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/bedtools"
transdecoder="/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/transdecoder"


for file in $(ls $small_gtf_dir); do
    base_name=$(basename "$file" | sed 's/\.[^.]*$//')
    print($base_name)
    # # Run bed tools (gtf to fasta)
    # bedtools getfasta -fi {input.genome} -bed config["scallop"]["output_file"] -fo config["bedtools"]["fasta_output"]

    # #transdecoder shit
    # gtf_to_alignment_gff3.pl config["scallop"]["output_file"] > {output.gff3}

    # TransDecoder.LongOrfs -t config["bedtools"]["fasta_output"] -O config["directories"]["transdecoder_dir"]
done