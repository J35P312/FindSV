import sys
import os
import argparse

#split the csq field
def fetch_genes(info_field):
    
    if"CSQ=" in info_field:
        csq_field= info_field.split("CSQ=")[-1]
    elif "ANN=" in info_field:
        csq_field= info_field.split("CSQ=")[-1]
    elif "EFF=" in info_field:
        csq_field= info_field.split("CSQ=")[-1]
    csq_field= csq_field.split(";")[0]
    genes=csq_field.split("|")
    for i in range(0,len(genes)):
        genes[i] = "|" + genes[i] + "|"
    return set(genes)    

parser = argparse.ArgumentParser("""the_Annotator: annotates gene specific information based on input txt,bed or tab files""")
parser.add_argument('--vcf',type=str,required=True,help="the path to the vcf file")
parser.add_argument('--folder',type=str,help="The path to the folder containing the files")
args, unknown = parser.parse_known_args()

gene_lists={}
for file in os.listdir(args.folder):
    if file.endswith(".txt") or file.endswith(".tab") or file.endswith(".bed"):
        gene_list_name=file.replace(".txt","").replace(".bed","").replace(".tab","")
        gene_list_name=gene_list_name.split("/")[-1]
        
        if not gene_list_name in gene_lists:
            gene_lists[gene_list_name] ={}
        for line in open(os.path.join(args.folder,file)):
            content=line.strip().split("\t")
            if len(content) > 1:
                gene_lists[gene_list_name]["|"+content[0]+"|"]=content[1]
            else:
                gene_lists[gene_list_name]["|"+content[0]+"|"]="None"

for line in open(args.vcf):
    if line[0] == "#" and line[1] == "#":
        print(line.strip())
    elif line[0] == "#":
        for gene_list in gene_lists:
          print  "##INFO=<ID={},Number=1,Type=String,Description=\"A gene list key\">".format(gene_list)
        print(line.strip())
    else:
        
        content=line.strip().split("\t")
        
        genes=fetch_genes(content[7])
        found_tags={}
        for gene in genes:
            for gene_list in gene_lists:
                if gene in gene_lists[gene_list]:
                    if gene_lists[gene_list][gene] == "None":
                        pass
                    else:
                        if not gene_list in found_tags:
                            found_tags[gene_list]="Found"
                             
                        found_tags[gene_list] += "|"+gene_lists[gene_list][gene]
        gene_tags=""
        for tag in found_tags:                
            gene_tags +=  tag+"=" + found_tags[tag]
        
        if len(found_tags) > 0:
            content[7] += ";" + gene_tags
        print "\t".join(content)
