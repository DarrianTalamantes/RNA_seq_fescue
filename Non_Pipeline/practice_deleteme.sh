# setting variables
genome="/scratch/drt83172/Wallace_lab/RNA_SEQ/Genome/tall_fescuev0.1.fa"
small_gtf_dir="/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/small_gtfs"
bedtools_dir="/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/bedtools"
gff3_dir="/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/gff3"
transdecoder="/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/transdecoder"


for file in $(ls $small_gtf_dir | grep "dupped" | grep "PxH_NegxPos"); do
    base_name=$(basename "$file" | sed 's/\.[^.]*$//')
    echo $base_name
    # # Run bed tools (gtf to fasta)
    bedtools getfasta -fi $genome -bed $small_gtf_dir/$base_name.gtf -fo $bedtools_dir/$base_name.fa

    # #transdecoder shit
    gtf_to_alignment_gff3.pl $small_gtf_dir/$base_name.gtf  > $gff3_dir/$base_name.gff3

    TransDecoder.LongOrfs -t $bedtools_dir/$base_name.fa -O $transdecoder
    
    TransDecoder.Predict -t $bedtools_dir/$base_name.fa -O $transdecoder
    # mv {params.pep_file_name} {params.output_dir}
    # mv {params.fasta_gff3_name} {params.output_dir}
    # mv {params.cds_name} {params.output_dir}
    # mv {params.bed_name} {params.output_dir}


done


