# Author: Darrian Talamantes
# Goal: For loop that downloads all the SRA reads that I want.
> download_commands.sh
for SRA in $(cat SRA_List.txt | cut -d "," -f 1) 
do
	echo 	"Processing $SRA" 
	fasterq-dump --threads 8 --progress --verbose $SRA -o $SRA 
	mv ${SRA}* /media/drt06/CTE_RNAseq/Extra_Seq
done

# parallel --verbose --jobs 1 < download_commands.sh 
