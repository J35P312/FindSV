import subprocess
import argparse

def main(args,parser):
    header={};
    header["ALT"]={}
    header["INFO"]={}
    header["FILTER"]={}
    header["FORMAT"]={}

    subheader={}
    columns=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"]
    
    print("##fileformat=VCFv4.1")
    print("##source=FindSV")
    for vcf in args.vcf:
        with open(vcf) as text_file:
            for line in text_file:
                if(line[0] == "#"):
                    if("#CHROM\tPOS" in line):
                        vcf_columns=line.strip().split("\t")
                        for column in vcf_columns:
                            if not column in columns:
                                columns.append(column.strip())
                        break
                    
                    elif line[0] == line[1] and line[0] == "#" and "=" in line:
                        if("=<ID=" in line):
                            field=line.split("=")[2].split(",")[0]
                            key= line.strip("#").split("=")[0]
                            if not key in header:
                                header[key]={field:line}
                            else:
                                header[key].update({field:line})
                        elif not "##source=" in line and not "##file" in line and not "##reference=" in line:
                            key= line.strip("#").split("=")[0]
                            if not key in subheader:
                                subheader.update({key:line})
                else:
                    break

    #print the mandatory header lines in the correct order
    for entry in sorted(header["ALT"]):
        print(header["ALT"][entry].strip())
    del header["ALT"]
    for entry in sorted(header["INFO"]):
        print(header["INFO"][entry].strip())
    del header["INFO"]
    for entry in sorted(header["FILTER"]):
        print(header["FILTER"][entry].strip())
    del header["FILTER"]

    for entry in sorted(header["FORMAT"]):
        print(header["FORMAT"][entry].strip())
    del header["FORMAT"]
    #print the other lines in lexiographic order
    for key in sorted(header):
        for entry in sorted(header[key]):
            print(header[key][entry].strip())


    #print subheaders
    for entry in sorted(subheader):
        print(subheader[entry].strip())

    print( "\t".join(columns) )


    for vcf in args.vcf:
        with open(vcf) as text_file:
            for line in text_file:
                if not (line[0] == "#"):
                    print(line.strip())

if __name__ == '__main__':
    parser = argparse.ArgumentParser("""Accepts multiple vcf files, merges them and prints them to stdout""")
    parser.add_argument('--vcf',nargs='+',type=str,required=True,help="the vcf files")
    args, unknown = parser.parse_known_args()
    main(args,parser)
