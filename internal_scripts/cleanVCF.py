import sys
import argparse

parser = argparse.ArgumentParser("""removes all vcf entries that are not marked pass or .""")
parser.add_argument('--vcf',type=str,required=True,help="the path to the vcf file")
args, unknown = parser.parse_known_args()

for line in open(args.vcf):
    
    if(line[0] == "#"):
        print(line.strip())
    else:
        content=line.split("\t")
        if content[6] == "PASS" or content[6] == "." or content[6] == "pass" or content[6] == "Pass":
            print(line.strip())

