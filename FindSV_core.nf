#!/usr/bin/env nextflow

if(params.folder){
    print "analysing all bam files in ${params.folder}\n"
    bam_files=Channel.fromPath("${params.folder}/*.bam")
    bai_files=Channel.fromPath("${params.folder}/*.bam.bai")
    if( bai_files.toList().length() == 0 ){
        bai_files=Channel.fromPath("${params.folder}/*.bai")
    }else{
        bai_files=Channel.fromPath("${params.folder}/*.bam.bai")
    }

}else if(params.bam){
    //first get the bam files, and check if all files exists   
    //bam_files= Channel.from( params.bam.splitCsv() )
    Channel.from( params.bam.splitCsv()).subscribe{
        if(!file(it).exists()) exit 1, "Missing bam:${it}, use either --bam to analyse a bam file, or --bam and--vcf to annotate a vcf file"
    }
    bam_files=Channel.from(params.bam.splitCsv()).map{
        line ->
        bam = file(line)
        [ bam ]
  }
    
    
    //then search for bai index, try both .bam.bai and .bai, stop if none is found
    Channel.from( params.bam.splitCsv() ).subscribe{
        if(file(it.replaceFirst(/.bam/,".bam.bai")).exists() ){
            println it.replaceFirst(/.bam/,".bam.bai")
        }else if(file(it.replaceFirst(/.bam/,".bai")).exists() ){
            println it.replaceFirst(/.bam/,".bai")
        }else{
            println "Missing bai:${it}, each bam must be indexed"
            exit 1
        }
    }
    
   //all bam file are indxed, now we can create a bam index channel
   bai_files=Channel.from(params.bam.splitCsv()).map {
        line ->
        if(file(line.replaceFirst(/.bam/,".bam.bai")).exists() ){
            bai = file(line.replaceFirst(/.bam/,".bam.bai"))    
        }else if(file(line.replaceFirst(/.bam/,".bai")).exists() ){
            bai = file(line.replaceFirst(/.bam/,".bai"))
        }
        [ bai ]
   }
   
}else{
    print "usage: nextflow FindSV_core.nf [--folder/--bam] --working_dir output_directory -c config_file\n"
    print "--bam STR,	analyse a bam file, the bam files is asumed to be indexed\nAnalyse multiple bam files by separating the path of each bam files by ,"
    print "--folder STR,	analyse all bam files in the given folder, the bam fils are assumed to be indexed\n"
    print "-c STR,	the config file generated using the setup.py script\n"
    print "--working_dir STR,	the output directory, here all the vcf files will end up\n"
    print "--vcf STR,	if used together with bam: the vcf file is annotated using information from the bam file\n if used with folder: not yet supported"
    exit 1
}


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

frequency_filter_exec= file("${params.frequency_filter_path}")

assemblatron_exec= file("${params.FindSV_home}/TIDDIT/variant_assembly_filter/assemblatron.nf")
assemblatron_conf_file = file("${params.FindSV_home}/TIDDIT/variant_assembly_filter/slurm.config")

TIDDIT_bam=Channel.create()
CNVnator_bam=Channel.create()
combine_bam=Channel.create()
annotation_bam=Channel.create()

Channel
       	.from bam_files
        .separate( TIDDIT_bam, CNVnator_bam, combine_bam,annotation_bam) { a -> [a, a, a,a,a] }


//perform variant calling if the input is a bam fle
if(!params.vcf){

    process TIDDIT {
        publishDir "${params.working_dir}", mode: 'copy', overwrite: true
        
        cpus 1
        
        input:
        file bam_file from TIDDIT_bam
    
        output:
        file "${bam_file.baseName}_FT_inter_chr_events.vcf" into TIDDIT_inter_vcf_files
        file "${bam_file.baseName}_FT_intra_chr_events.vcf" into TIDDIT_intra_vcf_files
    
        script:
        """
        ${TIDDIT_exec_file} --sv -b ${bam_file} -p ${params.TIDDIT_pairs} -q ${params.TIDDIT_q} -o ${bam_file.baseName}_FT
        """
    }

    process CNVnator {
        publishDir "${params.working_dir}", mode: 'copy', overwrite: true
        
        cpus 1

        input:
        file bam_file from CNVnator_bam
    
        output: 
            file  "${bam_file.baseName}_CNVnator.vcf" into CNVnator_vcf_files
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
        publishDir "${params.working_dir}", mode: 'copy', overwrite: true
        
        cpus 1

        input:
        file bam_file from combine_bam
        file CNVnator_vcf from CNVnator_vcf_files
        file TIDDIT_inter_vcf from TIDDIT_inter_vcf_files
        file TIDDIT_intra_vcf from TIDDIT_intra_vcf_files

        output: 
            file "${bam_file.baseName}_CombinedCalls.vcf" into combined_vcf_files
	    
	    
	    script:
	    
        """
        python ${SVDB_exec_file} --merge --no_var --pass_only --no_intra --overlap 0.7 --bnd_distance 2500 --vcf ${TIDDIT_inter_vcf} ${TIDDIT_intra_vcf} ${CNVnator_vcf} > merged.unsorted.vcf
        
        python ${contig_sort_exec_file} --vcf merged.unsorted.vcf --bam ${bam_file} > ${bam_file.baseName}_CombinedCalls.vcf
        rm merged.unsorted.vcf
        """
        
    }

    vcf_files=combined_vcf_files
    
}else{
	if(params.folder){
		vcf_files=Channel.fromPath("${params.vcf}")
		Channel.fromPath("${params.vcf}").subscribe{
            if(!file(it).exists()) exit 1, "Missing vcf:${it}, use either --bam to analyse a bam file, or --bam and--vcf to annotate a vcf file"
        }
	}else if(params.bam){
        vcf_files= Channel.from(params.vcf.splitCsv()).map {
            line ->
            vcf = file(line)
            [ vcf ]
        }
		
	}
	
}

process annotate{
    publishDir "${params.working_dir}", mode: 'copy', overwrite: true
    
    cpus 1
    
    input:
        file vcf_file from vcf_files
        file bam_file from annotation_bam
        file bai_file from bai_files
        file assemblatron_exec
        
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

    if [ "" != ${genmod_rank_model_file} ]
    then
        genmod score -c ${genmod_rank_model_file} ${bam_file.baseName}_FindSV.vcf  > ${vcf_file}.tmp
        mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    fi
	
    python ${frequency_filter_exec} ${bam_file.baseName}_FindSV.vcf ${params.SVDB_limit} > ${vcf_file}.tmp
    mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf

    if ["OFF" != ${params.assemblatron} ]
    then 
	${params.FindSV_home}/nextflow ${assemblatron_exec} --genome ${params.reference_fasta} --vcf ${bam_file.baseName}_FindSV.vcf --bam ${bam_file} --working_dir assemblatron_output -c ${assemblatron_conf_file}
        python ${cleanVCF_exec} --vcf assemblatron_output/assemblator.vcf > ${bam_file.baseName}_FindSV.vcf
    fi

    """

}
