# Author: Darrian Talamantes
# Goal: rename the downloaded SRA's to the important name


data_folder="/scratch/drt83172/Wallace_lab/RNA_SEQ/Data"
List="/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/Lists/SRA_List.txt"


while IFS=, read -r SRA NAME
do
    echo "Processing $data_folder/${SRA}_1.fastq with additional value $data_folder/${NAME}_1.fastq"

    mv $data_folder/${NAME}_1.fastq $data_folder/${NAME}_R1.fastq
    mv $data_folder/${NAME}_2.fastq $data_folder/${NAME}_R2.fastq
        
done < "$List"

