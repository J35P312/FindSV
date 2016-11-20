import sys
import os
#accepts the path to a vep annotated vcf as the first command, prints all structural variatn that reside within a gene
inputFile=sys.argv[1];
for line in open(inputFile):
    
    if(line[0] == "#"):
        print(line.strip())
    else:
        kg_frq=0
        exac_frq=0
        if ";1000GAF=" in line:
            content=line.split("\t")
            frq=content[7].split(";1000GAF=")[-1]
            kg_frq=frq.split(";")[0]
            
        if ";AF=" in line:
            content=line.split("\t")
            frq=content[7].split(";AF=")[-1]
            exac_frq=frq.split(";")[0]
        
        if float(exac_frq) <= 0.01 and float(kg_frq) <= 0.01:
            print(line.strip())

