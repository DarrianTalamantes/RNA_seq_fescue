# =================================================================================================
#     Trimming Single End Reads
# =================================================================================================

# I need to do input and output for both foward and reverse data. I made a list 
rule trim_galore_pe:
    input:
        fwd= datadir + "/{pairs}R1.fastq.gz",
        rev= datadir + "/{pairs}R2.fastq.gz"
    output:
        fasta_fwd=trimmed + "/{pairs}R1.fq.gz",
        report_fwd=trimmed + "/{pairs}R1_trimming_report.txt",
        fasta_rev=trimmed + "/{pairs}R2.fq.gz",
        report_rev=trimmed + "/{pairs}R2_trimming_report.txt"
    threads: 16
    params:
        extra="--illumina -q 20"
    wrapper:
        "v3.9.0/bio/trim_galore/pe"

