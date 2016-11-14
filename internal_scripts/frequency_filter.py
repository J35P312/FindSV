#first argument is the input vcf, second argument is the frequency limit, prints all variants lower or equal to the limit
import sys
for line in open(sys.argv[1]):
	if "#" == line[0]:
		print line.strip()
		continue

	if ";FRQ=" in line:
		info=line.split("\t")[7]
		frequency=info.split("FRQ=")[-1].split(";")[0]
		if float(frequency) <= float(sys.argv[2]):
			print line.strip()
	else:
		print line.strip()
