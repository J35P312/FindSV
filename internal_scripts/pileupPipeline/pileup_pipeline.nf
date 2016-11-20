
params.bam=""
bam_file=file(params.bam)
if(!file(params.bam).exists()) exit 1, "Missing bam"

params.ref="/home/jesper/reference/hg19.fa"
params.working_dir="pileup_output"
params.exac=""
params.kg=""
params.sweref=""
params.svvcf=""

VEP_exec_file="variant_effect_predictor.pl"
frequency_script="/home/jesperei/snpPipe/exac_annotation_sqlite.py"
frequency_filter="/home/jesperei/snpPipe/pop_freq_filter.py"
SVannotate_script="/home/jesperei/snpPipe/SVannotate.py"


process mpileup{
    publishDir params.working_dir , mode: 'copy', overwrite: true

    input:
    file bam_file
    
    output:
    file ("${bam_file.baseName}.snps.vcf") into SNP_vcf
    
    """
    samtools mpileup -uf ${params.ref} ${params.bam}  | bcftools view -Nvcg - > ${bam_file.baseName}.snps.vcf
    """

}

process annotate{
    publishDir params.working_dir , mode: 'copy', overwrite: true
    
    input:  
    file SNP_vcf
    file bam_file

    output:
    file "${bam_file.baseName}.snps.annotated.filtered.vcf" into vep_vcf
    file "${bam_file.baseName}.snps.annotated.frequency.vcf" into frequency_vcf
    file "${bam_file.baseName}.snps.annotated.frequency.filtered.vcf" into final_vcf
    

    """
    ${VEP_exec_file} --cache --force_overwrite  -i ${SNP_vcf}  -o ${bam_file.baseName}.snps.annotated.vcf --port 3337 --vcf --per_gene --format vcf -q
    grep -E "#|splice|HIGH|MODERATE|LOW" ${bam_file.baseName}.snps.annotated.vcf > ${bam_file.baseName}.snps.annotated.filtered.vcf

    cp ${bam_file.baseName}.snps.annotated.filtered.vcf ${bam_file.baseName}.snps.annotated.frequency.vcf

    if [ "" != ${params.exac} ]
    then
        cp ${bam_file.baseName}.snps.annotated.frequency.vcf ${bam_file.baseName}.snps.annotated.frequency.vcf.tmp
        python ${frequency_script} --vcf ${bam_file.baseName}.snps.annotated.frequency.vcf.tmp --exac ${params.exac} > ${bam_file.baseName}.snps.annotated.frequency.vcf
    fi

    if [ "" != ${params.kg} ]
    then
        cp ${bam_file.baseName}.snps.annotated.frequency.vcf ${bam_file.baseName}.snps.annotated.frequency.vcf.tmp
        python ${frequency_script} --vcf ${bam_file.baseName}.snps.annotated.frequency.vcf.tmp --exac ${params.kg} --tag  1000GAF > ${bam_file.baseName}.snps.annotated.frequency.vcf
    fi
    
    if [ "" != ${params.sweref} ]
    then
        cp ${bam_file.baseName}.snps.annotated.frequency.vcf ${bam_file.baseName}.snps.annotated.frequency.vcf.tmp
        python ${frequency_script} --vcf ${bam_file.baseName}.snps.annotated.frequency.vcf.tmp --exac ${params.sweref} > ${bam_file.baseName}.snps.annotated.frequency.vcf        
    fi 
    
    if [ "" != ${params.svvcf} ]
    then
        cp ${bam_file.baseName}.snps.annotated.frequency.vcf ${bam_file.baseName}.snps.annotated.frequency.vcf.tmp
        python ${SVannotate_script} --snpvcf ${bam_file.baseName}.snps.annotated.frequency.vcf.tmp --svvcf ${params.svvcf} > ${bam_file.baseName}.snps.annotated.frequency.vcf
    fi

    python ${frequency_filter} ${bam_file.baseName}.snps.annotated.frequency.vcf > ${bam_file.baseName}.snps.annotated.frequency.filtered.vcf    
    """    
}

