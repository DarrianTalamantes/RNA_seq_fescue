Steps to rerun Scripts.

1. First you have to  run the RNAseq.smk pipeline. It does not work on its own you must run each step seperatly but at least its easy to tell the order things are supposed to be ran and detects inputs and outputs.
2. You must run EnTAP with the results of the RNAseq.smk pipline. You use the Scallop2 output.
3. DeSeq2 can be run with the output of the RNAseq.smk pipeline. You use the Feature Counts output.