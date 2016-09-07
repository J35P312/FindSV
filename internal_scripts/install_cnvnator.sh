#!/bin/bash
cd $1
#install root and ad it to path

if [ `uname -o` == "GNU/Linux" ] ; 
then if [ `uname -p` == "x86_64" ] ; 
    then 
	gcc --version |grep 4.4 > /dev/null; 
	if [ $? -eq 0 ] ; 
	then
	    rootver=root_v5.34.34.Linux-slc6-x86_64-gcc4.4	   
	else
	    rootver=root_v5.34.34.Linux-ubuntu14-x86_64-gcc4.8
	fi 
    fi
fi

wget https://root.cern.ch/download/${rootver}.tar.gz
tar -xvf ${rootver}.tar.gz

cd root/bin/
source thisroot.sh
cd $1
#download and install cnvnator
wget https://github.com/abyzovlab/CNVnator/releases/download/v0.3.1/CNVnator_v0.3.1.zip
unzip CNVnator_v0.3.1.zip
cd CNVnator_v0.3.1
cd src/samtools
make 
cd ..
make
