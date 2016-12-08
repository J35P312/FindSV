#!/bin/sh
#script used for automated launcing of FindSV_core.nf.
#$1 bam_file
#$2 config_file
#$3 the output directory(aka publishdir)
#$4 the combined caller vcf file
#$5 tmp directory where the trace file is stored

echo  "SAMPLE_ID":$1
mkdir $5
./nextflow FindSV_core.nf --bam $1 -c $2 --working_dir $3 --vcf $4 -with-trace $5/trace.txt | tee $5/log.txt
