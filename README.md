FindSV
===========
FindSV is a structural variation pipeline written in nextflow and python. FindSV performs variant calling using TIDDIT and CNVnator, and Manta.
similar variants are merged, and the variants are annotated using VEP, frequency database, genmod, and custom annotation using the annotator script.

FindSV needs to be setup using the setup.sh script, this script will generate a config file aswell as a bash script for setting up the environment.

Run
===

    First initiate the FindSV environment:
        source ./FindSV_env.sh
        
    To analyse one bam file and put the output in the output_folder type:

        singularity exec FindSV.simg python FindSV.py --bam file.bam --output output_folder --config config_file

    To analyse a folder containing bam files type:

        singularity exec FindSV.simg python FindSV.py --folder input_folder --output output_folder --config config_file

Optionally, the pipeline may be run using the FindSV_core.nf script directly:
	
	The following command will analyse a bam file:
		singularity exec FindSV.simg nextflow FindSV_core.nf --bam file.bam --working_dir output -c config.conf
	An entire folder containing bam files could be analysed using this command
	    singularity exec FindSV.simg nextflow FindSV_core.nf --folder /the/bams/are/in/this/folder/ --working_dir output -c config.conf
	
		
Installation
============
Dependencies:

    Manta*
    variant effect predictor
    Singularity
    nextflow

*manta is optional, and will only be e run if activated

Install the dependencies, then download FindSV. 

    git clone https://github.com/J35P312/FindSV.git


next, you download the singularity image:

    cd FindSV
    singularity pull --name FindSV.simg shub://J35P312/FindSVSingularity

The image file  must  be stored in the FindSV folder!

lastly run the setup script:

    ./setup.sh

remember to download the vep cache files!
If you are using a server without internet connection, you will first have to download FindSV and pull the singularity image. 
Next you transfer these files to the server, and then you run the setup script.

Manta is turned of by default. To activate Manta, enter the config file, and set the  RunManta variable to anything except "FALSE". 
If you activate manta, you may have to increase the time limit. Most 30X bam files should be complete within 2 days of run time.

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
These files should contain two columns, the first column is the HGNC ID of the  gene, and the second column is a gene tag.
If a variant cover a tagged gene, that gene tag will be added to the vcf entry of that variant. As default, the OMIM id of any
affected gene will be added to the variant using the OMIM.txt tab file. Additionally, the OMIM Morbid genes will be flagged. These datassets were downloaded via the ENSEMBL biomart

Frequency database
==========
The frequency database is a vcf file. These vcf files are either multisample vcf files, or they contain OCC or FRQ tags as generated via SVDB. These files may be produced using SVDB, or any SV caller that performs multi-sample calling.
Additionally, the swegen SVDB files may be used as a database:

https://swefreq.nbis.se/

Genmod
========
Genmod is a tool for creating rankscores for variants. Visit the genmod website for more information.

https://github.com/moonso/genmod



