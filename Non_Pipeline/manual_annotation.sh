#!/bin/bash
#SBATCH -J small_gtf_ann
#SBATCH -p batch
#SBATCH --ntasks=8
#SBATCH --mem 40gb
#SBATCH -t 20:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/small_gtf_ann.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/small_gtf_ann.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu

echo "This JobID for this job is ${SLURM_JOB_ID}."
sleep 5
echo "Done."

sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU

ml TransDecoder/5.7.0-GCC-11.3.0
ml BEDTools/2.30.0-GCC-12.2.0
ml InterProScan/5.68-100.0-foss-2022a


# setting variables
genome="/scratch/drt83172/Wallace_lab/RNA_SEQ/Genome/tall_fescuev0.1.fa"
small_gtf_dir="/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/small_gtfs"
bedtools_dir="/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/bedtools"
gff3_dir="/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/gff3"
transdecoder="/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/transdecoder"
transdecoder_output="/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/transdecoder_output"
interpro="/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/interpro"
interpro_temp="/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/interpro_temp"


for file in $(ls $small_gtf_dir | grep "dupped"); do
    base_name=$(basename "$file" | sed 's/\.[^.]*$//')
    echo $base_name
    # # Run bed tools (gtf to fasta)
    bedtools getfasta -fi $genome -bed $small_gtf_dir/$base_name.gtf -fo $bedtools_dir/$base_name.fa

    # #transdecoder shit
    gtf_to_alignment_gff3.pl $small_gtf_dir/$base_name.gtf  > $gff3_dir/$base_name.gff3
    
    # Make directory for transdecoder to work
    mkdir $transdecoder/$base_name
    TransDecoder.LongOrfs -t $bedtools_dir/$base_name.fa -O $transdecoder/$base_name
    
    TransDecoder.Predict -t $bedtools_dir/$base_name.fa -O $transdecoder/$base_name
    mv $base_name.fa.transdecoder.pep $transdecoder_output
    mv $base_name.fa.transdecoder.gff3 $transdecoder_output
    mv $base_name.fa.transdecoder.cds $transdecoder_output
    mv $base_name.fa.transdecoder.bed $transdecoder_output
    
    # clean the pepfile
    sed 's/*//g' $transdecoder_output/$base_name.fa.transdecoder.pep > $transdecoder_output/$base_name.fa.transdecoder_clean.pep
    export JAVA_OPTS="-Xmx10G"
    interproscan.sh -cpu 8 -f TSV,GFF3 -goterms -b $interpro/$base_name -i $transdecoder_output/$base_name.fa.transdecoder_clean.pep -T $interpro_temp
done




