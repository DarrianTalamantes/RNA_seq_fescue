# This rule will run R feature counts.

rule feature_counts:
    input:
        gtf = config["scallop"]["output_file"],
        bams = expand(sep_bams + "{pairs}.bam", pairs=PAIRS)
    output:
        counts = config["directories"]["features"] + "feature_counts.txt"
    script:
        """
        feature_counts.R {input.gtf} {input.bams} {output.gtf}
        """