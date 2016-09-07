import sys
import argparse

parser = argparse.ArgumentParser("""stable sorting algorithm of genmod annotated vcf files""")
parser.add_argument('--vcf',type=str,required=True,help="the path to the vcf file")
args, unknown = parser.parse_known_args()

genmod_score_dictionary={}
for line in open(args.vcf):

    #print the header     
    if(line[0] == "#"):
        print(line.strip())
    else:

        #extract the rankcore as a float
        content=line.split("\t")
        RankScore=content[7].split("RankScore=")[-1];
        RankScore=RankScore.split(";")[0]
        RankScore=float(RankScore.split(":")[-1])
        #store in the container
        if not RankScore in genmod_score_dictionary:
            genmod_score_dictionary[RankScore]=[]
        genmod_score_dictionary[RankScore].append(line)

for score in sorted(genmod_score_dictionary,reverse=True):
    for line in genmod_score_dictionary[score]:
        print(line.strip())

