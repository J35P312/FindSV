import sys
import argparse
import subprocess
import os

def get_genes(line):
    genes=[]
    info=line.split("\t")[7]
    vep=info.split(";CSQ=")[-1].split(";")[0]
    for entry in vep.split(","):
        genes.append(entry.split("|")[3])
        
    return set(genes)
    
    
def main(args):
    SV_dictionary={}
    for line in open(args.svvcf):
        if line[0] == "#":
            continue
            
        genes=get_genes(line)
        
        var_id=line.split("\t")[2]
        for gene in genes:
            if not gene in SV_dictionary:
                SV_dictionary[gene]=[]
                
            SV_dictionary[gene].append(var_id)

    for line in open(args.snpvcf):
        if line[0] == "#":
            if not line[1] == "#":
                print "##INFO=<ID={},Number=1,Type=String,Description=\"variant id of compund sv\">".format(args.tag)
            print line.strip()
            continue
            
        genes=get_genes(line)
        output=line.split("\t")
        variants=[]
        for gene in genes:
            if gene in SV_dictionary:
                variants += SV_dictionary[gene]
        variants=list(set(variants))
        output[7] += ";{}={}".format(args.tag,"|".join(variants))
                
        print ( "\t".join(output).strip() )

parser = argparse.ArgumentParser("""this scripts annotates adds the id of structural variants to snps found in genes affected by SVs""")
parser.add_argument('--svvcf',type=str,required=True,help="the path to the sv vcf file")
parser.add_argument('--snpvcf',type=str,required=True,help="the path to the snp vcf file")
parser.add_argument('--tag',type=str,default="SVID",help="the vcf frequency tag(default = SVID)")
args, unknown = parser.parse_known_args()

main(args)

