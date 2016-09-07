#!/bin/bash
cd $1
#download and install cnvnator
wget https://github.com/abyzovlab/CNVnator/releases/download/v0.3.1/CNVnator_v0.3.1.zip
unzip CNVnator_v0.3.1.zip
cd CNVnator_v0.3.1
cd src/samtools
make 
cd ..
make
