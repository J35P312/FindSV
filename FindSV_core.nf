#!/usr/bin/env nextflow

params.bam="none"
bam_file = file(params.bam)


params.vcf="none"
vcf_file=file(params.vcf)

if(!bam_file.exists() && !vcf_file.exists()) exit 1, "Missing bam or vcf file. use either --bam to analyse a bam file, or --bam and--vcf to annotate a vcf file"

TIDDIT_exec_file = file( "${params.TIDDIT_path}" )
if(!TIDDIT_exec_file.exists()) exit 1, "Error: Missing TIDDIT executable, set the TIDDIT_path parameter in the config file"
SVDB_exec_file = file("${params.SVDB_script_path}")
if(!SVDB_exec_file.exists()) exit 1, "Error: Missing SVDB executable, set the SVDB_script_path parameter in the config file"

contig_sort_exec_file = file("${params.contig_sort_path}")
if(!contig_sort_exec_file.exists()) exit 1, "Error: Missing contig sort executable, set the contig sort parameter in the config file"

SVDB_file = file("${params.SVDB_path}")

genmod_rank_model_file = file("${params.genmod_rank_model_path}")

CNVnator_exec_file=params.CNVnator_path
//ROOT=file("${params.thisroot_path}")
CNVnator_reference_dir=file("${params.CNVnator_reference_dir_path}")

CNVnator2vcf=params.CNVnator2vcf_path
if("${params.CNVnator2vcf_path}" == "") exit 1, "Error: Missing cnvnator CNVnator2vcf script, set the CNVnator2vcf_path parameter in the config file"

VEP_exec_file=params.VEP_path

clear_vep_exec=file("${params.clear_vep_path}")
cleanVCF_exec=file("${params.cleanVCF_path}")

the_annotator_exec=file("${params.the_annotator_path}")
gene_keys_dir=file("${params.gene_keys_dir_path}")

//perform variant calling if the input is a bam fle
if(!vcf_file.exists()){

    process TIDDIT {
        publishDir "${params.working_dir}"
        
        cpus 1
        
        input:
        file bam_file
    
        output:
        file "${bam_file.baseName}_FT_inter_chr_events.vcf" into TIDDIT_inter_vcf
        file "${bam_file.baseName}_FT_intra_chr_events.vcf" into TIDDIT_intra_vcf
    
        script:
        """
        ${TIDDIT_exec_file} --sv -b ${bam_file} -p ${params.TIDDIT_pairs} -q ${params.TIDDIT_q} -o ${bam_file.baseName}_FT
        """
    }

    process CNVnator {
        publishDir "${params.working_dir}"
        
        cpus 1

        input:
        file bam_file
    
        output: 
            file  "${bam_file.baseName}_CNVnator.vcf" into CNVnator_vcf
        script:
        """
        
        ${CNVnator_exec_file} -root cnvnator.root -tree ${bam_file}
        ${CNVnator_exec_file} -root cnvnator.root -his ${params.CNVnator_bin_size} -d ${CNVnator_reference_dir}
        ${CNVnator_exec_file} -root cnvnator.root -stat ${params.CNVnator_bin_size} >> cnvnator.log
        ${CNVnator_exec_file} -root cnvnator.root -partition ${params.CNVnator_bin_size}
        ${CNVnator_exec_file} -root cnvnator.root -call ${params.CNVnator_bin_size} > ${bam_file.baseName}_CNVnator.out
        ${CNVnator2vcf} ${bam_file.baseName}_CNVnator.out >  ${bam_file.baseName}_CNVnator.vcf
        rm cnvnator.root
        """
    }

    process combine {
        publishDir "${params.working_dir}"
        
        cpus 1

        input:
        file bam_file
        file CNVnator_vcf
        file TIDDIT_inter_vcf
        file TIDDIT_intra_vcf

        output: 
            file "${bam_file.baseName}_CombinedCalls.vcf" into combined_vcf
	    
	    
	    script:
	    
        """
        python ${SVDB_exec_file} --merge --no_var --pass_only --no_intra --overlap 0.7 --bnd_distance 2500 --vcf ${TIDDIT_inter_vcf} ${TIDDIT_intra_vcf} ${CNVnator_vcf} > merged.unsorted.vcf
        
        python ${contig_sort_exec_file} --vcf merged.unsorted.vcf --bam ${bam_file} > ${bam_file.baseName}_CombinedCalls.vcf
        rm merged.unsorted.vcf
        """
        
    }

    vcf_file=combined_vcf
}



process annotate{
    publishDir "${params.working_dir}"
    
    cpus 1
    
    input:
        file vcf_file
        file bam_file
        
    output:
        file "${bam_file.baseName}_FindSV.vcf" into final_FindSV_vcf
        
    script:
    
    """
    ${VEP_exec_file} --cache --force_overwrite --poly b -i ${vcf_file}  -o ${vcf_file}.tmp --buffer_size 5 --port 3337 --vcf --per_gene --format vcf -q
    mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    python ${clear_vep_exec} ${bam_file.baseName}_FindSV.vcf > ${vcf_file}.tmp
    mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    
    python ${cleanVCF_exec} --vcf ${bam_file.baseName}_FindSV.vcf > ${vcf_file}.tmp
    mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    
    python ${SVDB_exec_file} --merge --overlap 1 --vcf --vcf ${bam_file.baseName}_FindSV.vcf > ${vcf_file}.tmp
    mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    python ${contig_sort_exec_file} --vcf ${bam_file.baseName}_FindSV.vcf --bam ${bam_file} > ${vcf_file}.tmp
    mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    
    
    if [ "" != ${gene_keys_dir} ] && [ "" != ${the_annotator_exec} ]
    then
        python ${the_annotator_exec} --folder ${gene_keys_dir} --vcf ${bam_file.baseName}_FindSV.vcf > ${vcf_file}.tmp
        mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    fi
    
    if [ "" != ${SVDB_file} ]
    then
        python ${SVDB_exec_file} --query --overlap ${params.SVDB_overlap} --bnd_distance ${params.SVDB_distance} --query_vcf ${bam_file.baseName}_FindSV.vcf --db ${SVDB_file} > ${vcf_file}.tmp
        mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    fi

    if [ "" != {genmod_rank_model_file} ]
    then
        genmod score -c ${genmod_rank_model_file} ${bam_file.baseName}_FindSV.vcf  > ${vcf_file}.tmp
        mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    fi

    """

}
