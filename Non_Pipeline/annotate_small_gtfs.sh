#!/bin/bash
#SBATCH -J small_gtf_ann
#SBATCH -p batch
#SBATCH --ntasks=24
#SBATCH --mem 120gb
#SBATCH -t 140:00:00
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

# This is to run after gtf_subsetter.py 


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

# TransDecoder.Predict -t {input.fasta} -O {params.output_dir}
# mv {params.pep_file_name} {params.output_dir}
# mv {params.fasta_gff3_name} {params.output_dir}
# mv {params.cds_name} {params.output_dir}
# mv {params.bed_name} {params.output_dir}

# # Cleans transdecoder output
# sed 's/*//g' {input.pep_file} > {output.pep_file_clean}

# # run interproscan
# export JAVA_OPTS="-Xmx50G"
# interproscan.sh -cpu 24 -f TSV,GFF3 -goterms -b /scratch/drt83172/Wallace_lab/RNA_SEQ/Annotation/interproscan_results -i /scratch/drt83172/Wallace_lab/RNA_SEQ/transcriptome/predicted_transcripts.fasta.transdecoder_dir/predicted_transcripts.fasta.clean.transdecoder.pep -T /scratch/drt83172/Wallace_lab/RNA_SEQ/Annotation/interproscan_temp







