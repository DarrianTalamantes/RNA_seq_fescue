# This rule will run R feature counts.
# Will need to change directoy to your scripts directy in execution

rule feature_counts:
    input:
        gtf = config["scallop"]["output_file"],
        bams = expand(config["directories"]["filtered_bams"] + "/{pairs}Aligned.sortedByCoord_filtered.out.bam", pairs=PAIRS)
    output:
        counts = config["directories"]["features"] + "feature_counts.txt"
    params:
        outpath = config["directories"]["features"]
    conda:
        "../Conda_Envs/R.yaml"
    # You can put R script path in config file so you dont have to change it here    
    shell:
        """
        if [ ! -d {params.outpath} ]; then 
        mkdir -p {params.outpath}; fi

        Rscript /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/Scripts/feature_counts.R {output.counts} {input.gtf} "{input.bams}" 
        """