#!/bin/sh
#script used for automated launcing of FindSV_core.nf.
#$1 bam_file
#$2 config_file
#$3 the output directory(aka publishdir)
#$4 the combined caller vcf file

echo  "SAMPLE_ID":$1
mkdir $3
FindSV_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "nextflow $FindSV_dir/FindSV_core.nf --bam $1 -c $2 --working_dir $3 --vcf $4 -with-trace $3/trace.txt"
nextflow $FindSV_dir/FindSV_core.nf --bam $1 -c $2 --working_dir $3 --vcf $4 -with-trace $3/trace.txt | tee $3/log.txt
