# =================================================================================================
#     Trimming Single End Reads
# =================================================================================================

# I need to do input and output for both foward and reverse data. I made a list 
rule trim_galore_pe:
    input:
        fastq = datadir + "/{sample}.fastq.gz"
    output:
        fasta = trimmed + "{sample}.fq.gz",
        report = trimmed + "{sample}.trimming_report.txt"
    threads: 1
    params:
        extra="--illumina -q 20",
    shell:
        """
        trim_galore --fastqc --gzip --output_dir {qcdir} {input.fastq}
        """

