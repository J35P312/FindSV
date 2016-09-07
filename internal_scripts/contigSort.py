import subprocess
import argparse

def main(args,parser):
    header=subprocess.check_output("samtools view -H {}".format(args.bam), shell=True)
    header=header.split("\n")
    chromosome_order=[]
    vcf_content={};
    meta_data="";
    for line in header:
        content=line.split("\t")
        if(content[0] == "@SQ"):
            chromosome=content[1].strip().split(":")[1]
            chromosome_order.append(chromosome)
            vcf_content[chromosome]=[];
    with open(args.vcf) as text_file:
        for line in text_file:
            if(line[0] == "#"):
                meta_data+=line
            else:
                chromosome=line.split("\t")[0];
                #if the caller output is correct according to the bam reference
                if chromosome in vcf_content:
                    vcf_content[chromosome].append(line);
                #add chr prefix if needed
                elif "chr"+chromosome in vcf_content:
                    vcf_content["chr"+chromosome].append("chr"+line);
                #remove chr if needed
                elif chromosome[3:] in vcf_content:
                    vcf_content[chromosome[3:]].append(line[3:]);
    print(meta_data.strip())
    for chromosome in chromosome_order:
        positions={}
        for line in vcf_content[chromosome]:
            pos=int(line.split("\t")[1])
            if pos in positions:
                positions[pos].append(line)
            else:
                positions[pos]=[line]
        for pos in sorted(positions):
            for variant in positions[pos]:
                print(variant.strip())

if __name__ == '__main__':
    parser = argparse.ArgumentParser("""Uses samtools to read the header data of a bam file, accepts a vcf and sorts the contigs of the vcf according to the contig order, the output is printed to stdout""")
    parser.add_argument('--vcf',type=str,required=True,help="the path to the vcf file")
    parser.add_argument('--bam',type=str,required=True,help="the path to the bam file")
    args, unknown = parser.parse_known_args()
    main(args,parser)
