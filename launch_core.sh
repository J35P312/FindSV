#!/bin/sh
#script used for automated launcing of FindSV_core.nf.
#$1 bam_file
#$2 config_file
#$3 the output directory(aka publishdir)

echo  "SAMPLE_ID":$1
./nextflow FindSV_core.nf --bam $1 -c $2 --working_dir $3 -with-trace $4/trace.txt | tee $4/log.txt

