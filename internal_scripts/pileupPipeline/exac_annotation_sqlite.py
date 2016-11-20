import sys
import argparse
import subprocess
import os
import fnmatch
import sqlite3


def main(args):

    conn = sqlite3.connect(':memory:')
    c = conn.cursor()
    A="CREATE TABLE SVDB (chr TEXT, pos INT,alt TEXT, frequency FLOAT)"
    c.execute(A)

    exac_db=[]
    for line in open(args.exac):
        if line[0] == "#":
            continue
        content=line.strip().split()
        if not "," in content[4]:
            exac_db.append([ content[0], content[1] , content[3] ,content[4] ])
        else:
            alt_list=content[3].split(",")
            frq_list=content[4].split(",")
            for i in range(0,len(alt_list)):
                exac_db.append([ content[0], content[1] , alt_list[i] , frq_list[i] ])
                
        if len(exac_db) > 100000:
            c.executemany('INSERT INTO SVDB VALUES (?,?,?,?)',exac_db)          
            exac_db=[]
    if exac_db:
        c.executemany('INSERT INTO SVDB VALUES (?,?,?,?)',exac_db)
                    
    del exac_db
    
    A="CREATE INDEX SNP ON SVDB (chr, pos, alt)"
    c.execute(A)
    conn.commit()
    
    for line in open(args.vcf):
        if line[0] == "#" and line[1] == "#":
            print line.strip()
        elif line[0] == "#":
            print "##INFO=<ID={},Number=1,Type=Float,Description=\"Estimated allele frequency in the range (0,1).\">".format(args.tag)
            print line.strip()
        else:
            FRQ=0
            annotation=";{}={}"
            content=line.split("\t")

            chromosome=content[0]
            position=content[1]
            end=position
            alt=content[4]
            #cadd annotation
            
            A='SELECT frequency FROM SVDB WHERE chr == \'{}\' AND pos == {} AND alt == \'{}\' '.format(chromosome,position,alt)
            
            hits=[]
            for hit in c.execute(A):
                FRQ = hit[0]
            
            
            

            content[7] += annotation.format(args.tag,FRQ)
            print("\t".join(content).strip())
	   		#popfreq
        
    conn.close()
parser = argparse.ArgumentParser("""this scripts annotates a SNP using CADD and popfreq, the popfreq and CADD file must be tabix indexed and tabbix must be installed""")
parser.add_argument('--vcf',type=str,required=True,help="the path to the vcf file")
parser.add_argument('--folder',type=str,help="used instead of vcf to annotate each vcf in a folder")
parser.add_argument('--exac',type=str,help="the path to the exac DB")
parser.add_argument('--tag',type=str,default="AF",help="the vcf frequency tag(default = AF)")
args, unknown = parser.parse_known_args()

if args.vcf:
    main(args)
elif args.folder:
    for root, dirnames, filenames in os.walk(args.folder):
            for filename in fnmatch.filter(filenames, '*.vcf'):
                bam_file=os.path.join(root, filename)
                args.vcf=bam_file
                main(args)
else:
    print("|>L3453 5|>3(1/=Y --vcf || --folder")
