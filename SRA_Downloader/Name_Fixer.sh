# Author: Darrian Talamantes
# Goal: rename the downloaded SRA's to the important name


data_folder="/scratch/drt83172/Wallace_lab/RNA_SEQ/Data"
List="/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/Lists/Fixed_Names_List.csv"


while IFS=, read -r SRA NAME
do
    echo "Processing $data_folder/${SRA} with additional value $data_folder/${NAME}"

    mv $data_folder/${SRA} $data_folder/${NAME}
    mv $data_folder/${SRA} $data_folder/${NAME}
        
done < "$List"

