FindSV
===========
FindSV is a structural variation pipeline written in nextflow and python. FindSV performs variant calling using TIDDIT and CNVnator.
similar variants are merged, and the variants are annotated using VEP, frequency database, genmod, and custom annotation using the annotator script.

FindSV needs to be setup using the setup.py script, this script will generate a config file aswell as a bash script for setting up the environment.

Run
===

    First initiate the FindSV environment:
        source ./FindSV_env.sh
        
    To analyse one bam file and put the output in the output_folder type:

        python FindSV.py --bam file.bam --output output_folder --config config_file

    To analyse a folder containing bam files type:

        python FindSV.py --folder input_folder --output output_folder --config config_file

	Optionally, the pipeline may be run using the FindSV_core.nf script directly:
		./nextflow FindSV_core.nf --bam file.bam --working_dir output -c config.conf


Installation
============
Dependencies:
    Conda
    cnvnator
    variant effect predictor
    bwa
    
After installing the dependencies, run the setup script:
    python setup.py
    
this script will ask a couple on questions, such as cnvnator path and reference directory path. answer all these questions to finnish the setup.
The cnvnator install scripts in the internal_scripts folder may be used to intall cnvnator.

After running the setup script, two files will be generated, the config file, and the FindSV_env.sh script.

Restart module
============
for info on how to restart samples analysed by the FindSV- pipeline, type:

        python FindSV.py --restart

Config file
=========
The config file is generated using the setup.py script. Generate the file and read the comments within 
this file for more information on the individual variables.

The config file template is found in the folder config_template. This template could be turned into a config file through manual configuration.

Gene_Keys
==========

add tab files to this folder to annotate genes.
These files should contain two columns, the first column is the gene, and the second column is a gene tag.
If a variant cover a tagged gene, that gene tag will be added to the vcf entry of that variant. As default, the OMIM id of any
affected gene will be added to the variant using the OMIM.txt tab file.

Frequency database
==========

Genmod
========


