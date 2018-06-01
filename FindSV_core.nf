#!/usr/bin/env nextflow

params.folder = ""
params.vcf=""

if(params.folder){
    print "analysing all bam files in ${params.folder}\n"
	
    bam_files=Channel.fromPath("${params.folder}/*.bam").map{
        line -> 
        ["${file(line).baseName}",file(line)]
	}

    bai_files=Channel.fromPath("${params.folder}/*.bai").map{
        line ->
        println "${file(line).baseName}"
        [ "${file(line).baseName}" , file(line) ]
    }
    
}else if(params.bam){
    //first get the bam files, and check if all files exists   
    //bam_files= Channel.from( params.bam.splitCsv() )
    Channel.from( params.bam.splitCsv()).subscribe{
        if(!file(it).exists()) exit 1, "Missing bam:${it}, use either --bam to analyse a bam file, or --bam and--vcf to annotate a vcf file"
    }
	bam_files=Channel.from(params.bam.splitCsv()).map{
        line ->
        ["${file(line).baseName}", file(line) ]
	}	
    

    //search for bai index, try both .bam.bai and .bai, stop if none is found
    bai_files=Channel.from(params.bam.splitCsv()).map{
        line ->
        if(file(line+".bai").exists() ){
            [ "${file(line).baseName}",file(line+".bai") ]   
        }else if(file(line.replaceFirst(/.bam/,".bai")).exists() ){
            [ "${file(line).baseName}", file(line.replaceFirst(/.bam/,".bai")) ]
        }else{
            println "Missing bai:${line}, each bam must be indexed"
            exit 1
        }
    }
    if (params.vcf){
        vcf_files=Channel.from(params.bam.splitCsv()).map{
            line ->
            if(file(line.replaceFirst(/.bam/,".bam.bai")).exists() ){
                bai = file(line.replaceFirst(/.bam/,".bam.bai"))
                [ file(line), file(line.replaceFirst(/.bam/,".bam.bai")) , file("${params.vcf}/${file(line).baseName}_CombinedCalls.vcf") ]   
            }else if(file(line.replaceFirst(/.bam/,".bai")).exists() ){
                [ file(line), file(line.replaceFirst(/.bam/,".bai")) , file("${params.vcf}/${file(line).baseName}_CombinedCalls.vcf") ]
            }
        } 
        
        Channel.from(params.bam.splitCsv()).map{
            line ->
            if(file(line.replaceFirst(/.bam/,".bam.bai")).exists() ){
                bai = file(line.replaceFirst(/.bam/,".bam.bai"))
                [ file(line), file(line.replaceFirst(/.bam/,".bam.bai")) , file("${params.vcf}/${file(line).baseName}_CombinedCalls.vcf") ]   
            }else if(file(line.replaceFirst(/.bam/,".bai")).exists() ){
                [ file(line), file(line.replaceFirst(/.bam/,".bai")) , file("${params.vcf}/${file(line).baseName}_CombinedCalls.vcf") ]
            }
        }.subscribe{println it}   
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

sample_data = bam_files.cross(bai_files).map{
      it ->  [it[0][0],it[0][1],it[1][1]] }

contig_sort_exec_file = file("${params.contig_sort_path}")
if(!contig_sort_exec_file.exists()) exit 1, "Error: Missing contig sort executable, set the contig sort parameter in the config file"

SVDB_file = file("${params.SVDB_path}")
pathogenic_db_file=file("${params.pathogenic_db_path}")
benign_db_file=file("${params.benign_db_path}")

genmod_rank_model_file = params.genmod_rank_model_path

CNVnator_reference_dir=file("${params.CNVnator_reference_dir_path}")

VEP_exec_file=params.VEP_path

clear_vep_exec=file("${params.clear_vep_path}")
cleanVCF_exec=file("${params.cleanVCF_path}")

the_annotator_exec=file("${params.the_annotator_path}")
gene_keys_dir=file("${params.gene_keys_dir_path}")

frequency_filter_exec= file("${params.frequency_filter_path}")

TIDDIT_bam=Channel.create()
CNVnator_bam=Channel.create()
manta_bam=Channel.create()
annotation_bam=Channel.create()
combined_bam=Channel.create()

Channel
       	.from sample_data
        .separate( TIDDIT_bam,annotation_bam,combined_bam, CNVnator_bam,manta_bam) { a -> [a, a, a, a,a] }
        
//perform variant calling if the input is a bam fle
if(!params.vcf){
    if ("${params.RunManta}" != "FALSE"){
        //Then assemble the variants
        process Manta {
            publishDir "${params.working_dir}", mode: 'copy', overwrite: true
            errorStrategy 'ignore'      
            tag { bam_file }
            scratch true
            cpus 2
    
            input:
            set ID,  file(bam_file), file(bai_file) from manta_bam
    
            output:
            set ID, "${bam_file.baseName}.diploidSV.vcf"  into manta_output
    
            script:
            """
            ${params.configManta} --normalBam ${bam_file} --reference ${params.genome} --runDir MANTA_DIR
            python MANTA_DIR/runWorkflow.py -m local -j 2
            gunzip -c MANTA_DIR/results/variants/diploidSV.vcf.gz  > ${bam_file.baseName}.diploidSV.vcf
            sed -ie 's/DUP:TANDEM/TDUP/g' ${bam_file.baseName}.diploidSV.vcf
            """
        }   
    }

    process TIDDIT {
        publishDir "${params.working_dir}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'      
        tag { bam_file }
        scratch true
    
        cpus 1
        
        input:
        set ID,  file(bam_file), file(bai_file) from TIDDIT_bam
    
        output:
        set ID, "${bam_file.baseName}.vcf" into TIDDIT_output
    
        script:
        """
        singularity exec ${params.FindSV_home}/FindSV.simg python /opt/TIDDIT/TIDDIT.py --sv --bam ${bam_file} -p ${params.TIDDIT_pairs} -q ${params.TIDDIT_q} -o ${bam_file.baseName} --ref ${params.genome}
        rm *.tab
        """
    }

    process CNVnator {
        publishDir "${params.working_dir}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'
        tag { bam_file }       
        scratch true

        cpus 1

        input:
        set ID,  file(bam_file), file(bai_file) from CNVnator_bam
    
        output: 
        set ID, "${bam_file.baseName}_CNVnator.vcf" into CNVnator_output
            
        script:
        """
        
        singularity exec ${params.FindSV_home}/FindSV.simg /opt/CNVnator_v0.3.3/src/cnvnator -root cnvnator.root -tree ${bam_file}
        singularity exec ${params.FindSV_home}/FindSV.simg /opt/CNVnator_v0.3.3/src/cnvnator -root cnvnator.root -his ${params.CNVnator_bin_size} -d ${CNVnator_reference_dir}
        singularity exec ${params.FindSV_home}/FindSV.simg /opt/CNVnator_v0.3.3/src/cnvnator -root cnvnator.root -stat ${params.CNVnator_bin_size} >> cnvnator.log
        singularity exec ${params.FindSV_home}/FindSV.simg /opt/CNVnator_v0.3.3/src/cnvnator -root cnvnator.root -partition ${params.CNVnator_bin_size}
        singularity exec ${params.FindSV_home}/FindSV.simg /opt/CNVnator_v0.3.3/src/cnvnator -root cnvnator.root -call ${params.CNVnator_bin_size} > ${bam_file.baseName}_CNVnator.out
        singularity exec ${params.FindSV_home}/FindSV.simg /opt/CNVnator_v0.3.3/cnvnator2VCF.pl ${bam_file.baseName}_CNVnator.out >  ${bam_file.baseName}_CNVnator.vcf
        rm cnvnator.root
        """
    }
    


    if ("${params.RunManta}" != "FALSE"){

        combined_TIDDIT_CNVnator = TIDDIT_output.cross(CNVnator_output).map{
            it ->  [it[0][0],it[0][1],it[1][1]]
        }

        combined_data = combined_TIDDIT_CNVnator.cross(manta_output).map{
            it ->  [it[0][0],it[0][1],it[0][2],it[1][1]]
        }
    }else{
        combined_data = TIDDIT_output.cross(CNVnator_output).map{
            it ->  [it[0][0],it[0][1],it[1][1],file(params.genome)]
        }

    }    

    combined_bam = combined_data.cross(combined_bam).map{
        it ->  [it[0][0],it[1][1],it[1][2],it[0][1],it[0][2],it[0][3] ]
    }

    process combine {
        publishDir "${params.working_dir}", mode: 'copy', overwrite: true
        errorStrategy 'ignore'        
        tag { bam_file }
        cpus 1

        input:
        set ID,file(bam_file), file(bai_file), file(TIDDIT_vcf), file(CNVnator_vcf), file(manta_vcf) from combined_bam

        output: 
        set ID, file("${bam_file.baseName}_CombinedCalls.vcf") into combined_FindSV
	    
	    
	    script:

        if ("${params.RunManta}" != "FALSE"){
            """
            singularity exec ${params.FindSV_home}/FindSV.simg svdb --merge --no_var --no_intra --same_order --overlap 0.7 --bnd_distance 2500 --vcf ${TIDDIT_vcf} ${manta_vcf} > PE_signals.vcf
            singularity exec ${params.FindSV_home}/FindSV.simg svdb --merge --pass_only --no_intra --same_order --overlap 0.7 --bnd_distance 2500 --vcf PE_signals.vcf ${CNVnator_vcf} > merged.unsorted.vcf
            singularity exec ${params.FindSV_home}/FindSV.simg python ${contig_sort_exec_file} --vcf merged.unsorted.vcf --bam ${bam_file} > ${bam_file.baseName}_CombinedCalls.vcf
            rm merged.unsorted.vcf
            """
        }else{
            """
            singularity exec ${params.FindSV_home}/FindSV.simg svdb --merge --pass_only --same_order --no_intra --overlap 0.7 --bnd_distance 2500 --vcf ${TIDDIT_vcf} ${CNVnator_vcf} > merged.unsorted.vcf
            singularity exec ${params.FindSV_home}/FindSV.simg python ${contig_sort_exec_file} --vcf merged.unsorted.vcf --bam ${bam_file} > ${bam_file.baseName}_CombinedCalls.vcf
            rm merged.unsorted.vcf
            """
        }


    }

    vcf_files=annotation_bam.cross(combined_FindSV).map{
        it ->  [it[0][1],it[0][2],it[1][1]]
    }
    
}

process annotate{
    publishDir "${params.working_dir}", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    tag { bam_file } 

    cpus 1
    
    input:
    set file(bam_file), file(bai_file), file(vcf_file) from vcf_files

    output:
        file "${bam_file.baseName}_FindSV.vcf" into final_FindSV_vcf
        
    script:
    
    """
    ${VEP_exec_file} -i ${vcf_file}  -o ${vcf_file}.tmp ${params.vep_args}
    mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    singularity exec ${params.FindSV_home}/FindSV.simg python ${clear_vep_exec} ${bam_file.baseName}_FindSV.vcf > ${vcf_file}.tmp
    mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    
    singularity exec ${params.FindSV_home}/FindSV.simg python ${cleanVCF_exec} --vcf ${bam_file.baseName}_FindSV.vcf > ${vcf_file}.tmp
    mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    
    singularity exec ${params.FindSV_home}/FindSV.simg svdb --merge --overlap 0.9 --vcf ${bam_file.baseName}_FindSV.vcf > ${vcf_file}.tmp
    mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    singularity exec ${params.FindSV_home}/FindSV.simg python ${contig_sort_exec_file} --vcf ${bam_file.baseName}_FindSV.vcf --bam ${bam_file} > ${vcf_file}.tmp
    mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    
    
    if [ "" != ${gene_keys_dir} ] && [ "" != ${the_annotator_exec} ]
    then
        singularity exec ${params.FindSV_home}/FindSV.simg python ${the_annotator_exec} --folder ${gene_keys_dir} --vcf ${bam_file.baseName}_FindSV.vcf > ${vcf_file}.tmp
        mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    fi
    
    if [ "" != ${params.SVDB_path} ]
    then
        singularity exec ${params.FindSV_home}/FindSV.simg svdb --query --overlap ${params.SVDB_overlap} --bnd_distance ${params.SVDB_distance} --query_vcf ${bam_file.baseName}_FindSV.vcf --db ${SVDB_file} > ${vcf_file}.tmp
        mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    fi

    if [ "" != ${params.pathogenic_db_path} ]
    then 
        singularity exec ${params.FindSV_home}/FindSV.simg svdb --query --overlap ${params.pathogenic_db_overlap} --bnd_distance ${params.pathogenic_db_distance} --query_vcf ${bam_file.baseName}_FindSV.vcf --db ${pathogenic_db_file} --hit_tag PATHOGENIC > ${vcf_file}.tmp
        mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    fi


    if [ "" != ${genmod_rank_model_file} ]
    then
        singularity exec ${params.FindSV_home}/FindSV.simg genmod score -c ${genmod_rank_model_file} ${bam_file.baseName}_FindSV.vcf  > ${vcf_file}.tmp
        mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf
    fi
	
    singularity exec ${params.FindSV_home}/FindSV.simg python ${frequency_filter_exec} ${bam_file.baseName}_FindSV.vcf ${params.SVDB_limit} > ${vcf_file}.tmp
    mv ${vcf_file}.tmp ${bam_file.baseName}_FindSV.vcf

    if [ "" != ${params.benign_db_path} ]
    then 
        singularity exec ${params.FindSV_home}/FindSV.simg svdb --query --overlap ${params.benign_db_overlap} --bnd_distance ${params.benign_db_distance} --query_vcf ${bam_file.baseName}_FindSV.vcf --db ${benign_db_file} --hit_tag BENIGN > ${vcf_file}.tmp
        grep -v ";BENIGN=1" ${vcf_file}.tmp > ${bam_file.baseName}_FindSV.vcf
        rm ${vcf_file}.tmp
    fi

    """

}
