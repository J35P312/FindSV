import sys
# removes duplicate gene entries so that there will only be one VEP entry per gene within each variant
#accepts one vcf file as command line argument
for line in open(sys.argv[1]):
    if "#" == line[0]:
        print line.strip()
        continue

    if not "CSQ=" in line or "intergenic" in line:
        print line.strip()
        continue

    CSQ= line.strip().split("\t")[7].split("CSQ=")[-1].split(";")[0]
    vcf_line=line.split("\t")
    info_field=vcf_line[7].strip().split(";")
 

    CSQ_fields=CSQ.split(",")
    genes={}
    for field in CSQ_fields:
        content= field.split("|")
        if not content[3] in genes:
            genes[content[3]] = field
        else:

            if "|HIGH|" in genes[content[3]]:
                pass
            elif "|MODIFIER|" in genes[content[3]]:
                genes[content[3]] = field
            elif "|MODERATE|" in genes[content[3]] and "|HIGH|" in field:
                genes[content[3]] = field
            elif "|MODIFIER|" in genes[content[3]] or "|LOW|" in genes[content[3]] and "|MODERATE|" in field:
                genes[content[3]] = field
    field_list=[]
    for gene in genes:
        field_list.append(genes[gene])
    csq="CSQ=" + ",".join(field_list)
    i=0
    while i < len(info_field):
        if info_field[i].startswith("CSQ="):
            info_field[i] =  csq

        i +=1     
    vcf_line[7]=";".join(info_field)
    print "\t".join(vcf_line).strip()
