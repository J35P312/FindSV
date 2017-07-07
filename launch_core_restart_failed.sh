#!/bin/sh
#script used for automated launcing of FindSV_core.nf.
#sample["bam_call"],sample["bam_annotate"],args.config,args.output,sample["vcf"],temp_dir

FindSV_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo  "SAMPLE_ID":$2
if [ "NONE" != "$2" ] && [ "NONE" != "$1"]
    then
        $FindSV_dir/nextflow $FindSV_dir/FindSV_core.nf --bam $2 -c $3 --working_dir $4 --vcf $5 -with-trace $6/trace_annotate.txt | tee $6/log_annotation.txt &
    elif [ "NONE" != "$2" ]
        then
            $FindSV_dir/nextflow $FindSV_dir/FindSV_core.nf --bam $2 -c $3 --working_dir $4 --vcf $5 -with-trace $6/trace_annotate.txt | tee $6/log_annotation.txt
    else
       > $6/log_annotation.txt
fi
echo  "SAMPLE_ID":$1

if [ "NONE" != "$1" ]
    then
    $FindSV_dir/nextflow $FindSV_dir/FindSV_core.nf --bam $1 -c $3 --working_dir $4 -with-trace $6/trace_call.txt | tee $6/log_calling.txt
    else
        > $6/log_calling.txt
fi


cat $6/log_calling.txt $6/log_annotation.txt > $6/log.txt
cat $6/trace_call.txt $6/trace_annotate.txt > $6/trace.txt

