import os
import subprocess

programDirectory = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(programDirectory,"config_template/FindSV.config"), 'r') as myfile:
    template=myfile.read()

valid_response=["1","2"]
while True:
    print "print 1 for local executor, or 2 for slurm"
    selection=raw_input()
    if not selection in valid_response:
        print "invalid option, print 1 for local executor, or 2 for slurm"
    else:

        if selection == "1":
            selection = "local"
        else:   
            selection = "slurm"
        break
        
template=template.replace("{executor}", "\'{}\'".format(selection) )
if selection == "slurm":
    print "enter slurm account, or press enter"
    selection=raw_input()
    
template=template.replace("{account}", "{}".format(selection) )

print "enter the standard output directory"
selection=raw_input()
template=template.replace("{working_dir}", "\'{}\'".format(selection) )

print "installing and setting up TIDDIT"
command=["{} {}".format(os.path.join(programDirectory,"internal_scripts/install_FT.sh"),programDirectory)]
tmp=subprocess.check_output(command,shell = True)
template=template.replace("{TIDDIT_path}", "\'{}\'".format(os.path.join(programDirectory,"TIDDIT/bin/TIDDIT")) )

print "Setting up manta"
print "Set the manta configManta path, the path is set to configManta if left blank"
selection=raw_input()
if selection == "":
    selection = "configManta"
template=template.replace("configManta}", "\'{}\'".format(selection) )

print "Set the path to the reference genome(needs to be indexed using bwa and samtools 0.19)"
selection=raw_input()
template=template.replace("{genome}", "\'{}\'".format(selection) )

print "add cnvnator path, the path is set to cnvnator if left blank"
selection=raw_input()
if selection == "":
    selection = "cnvnator"
template=template.replace("{CNVnator_path}", "\'{}\'".format(selection) )
  
print "if the ROOTSYS variable is not set, add the path to the thisroot.sh script inside the root bin folder, leave blank otherwise(open another terminal and print echo $ROOTSYS to check)"
selection=raw_input()
template=template.replace("{thisroot_path}", "\'{}\'".format(selection) )
thisroot=selection

print "add the path of the directory containing reference files"
selection=raw_input()
template=template.replace("{CNVnator_reference_dir_path}", "\'{}\'".format(selection) )

print "add the cnvnator2VCF.pl script path, the path is set to cnvnator2VCF.pl if left blank"
selection=raw_input()
if selection == "":
    selection = "cnvnator2VCF.pl"
template=template.replace("{CNVnator2vcf_path}", "\'{}\'".format(selection) )


print "setting up internal pipeline scripts"
template=template.replace("{contig_sort_path}", "\'{}\'".format( os.path.join(programDirectory,"internal_scripts/contigSort.py") ) )
template=template.replace("{clear_vep_path}", "\'{}\'".format( os.path.join(programDirectory,"internal_scripts/clear_vep.py") ) )
template=template.replace("{cleanVCF_path}", "\'{}\'".format( os.path.join(programDirectory,"internal_scripts/cleanVCF.py") ) )
template=template.replace("{the_annotator_path}", "\'{}\'".format( os.path.join(programDirectory,"internal_scripts/the_annotator.py") ) )
template=template.replace("{gene_keys_dir_path}", "\'{}\'".format( os.path.join(programDirectory,"gene_keys") ) )
template=template.replace("{frequency_filter_path}", "\'{}\'".format( os.path.join(programDirectory,"internal_scripts/frequency_filter.py") ) )
template=template.replace("{FindSV_home_path}", "\'{}\'".format( os.path.join(programDirectory) ) )
template=template.replace("{run_fermikit_path}", "\'{}\'".format( os.path.join(programDirectory,"internal_scripts/run_fermikit.py") ) )

print "add the variant_effect_predictor.pl script path, the path is set to variant_effect_predictor.pl if left blank"
print "remember to download the VEP cache file! more info is found on the VEP ENSMBLE website"
selection=raw_input()
if selection == "":
    selection = "variant_effect_predictor.pl"
template=template.replace("{VEP_path}", "\'{}\'".format(selection) )

print "intalling SVDB"
template=template.replace("{SVDB_script_path}", "\'{}\'".format( os.path.join(programDirectory,"SVDB/SVDB.py") ) )
os.system("git clone https://github.com/SciLifeLab/SVDB.git")

print "add the path of an SVDB sqlite database file"
selection=raw_input()
template=template.replace("{SVDB_path}", "\'{}\'".format(selection) )

print "add the path of vcf database file containing known pathogenic variants(or press enter to continue)"
selection=raw_input()
if selection == "":
    selection = "\"\""
template=template.replace("{pathogenic_db_path}", "\'{}\'".format(selection) )


print "add the path of vcf database file containing known benign variants(or press enter to continue)"
selection=raw_input()
if selection == "":
    selection = "\"\""
template=template.replace("{benign_db_path}", "\'{}\'".format(selection) )

print "add the genmod ini file path, if left blank, the default genmod_SV.txt file will be used"
selection=raw_input()
if selection == "":
    selection = os.path.join(programDirectory,"genmod_SV.txt")
template=template.replace("{genmod_rank_model_path}", "\'{}\'".format(selection) )

print "enter the filename of the new config file(or leave blank to set the config file to FindSV.conf)"
selection=raw_input()
if selection == "":
    selection = "FindSV.conf"
f= open(selection, "w")
f.write(template)
f.close()

print "creating FindSV environment script"
print "modules: print uppmax if you are using uppmax, print a line of each module to use, or leave empty to skip modules"
print "example: bioinfo-tools samtools CNVnator vep, to load the modules bioinfo-tools, sammtools, CNVnator and vep"
selection=raw_input()
if selection == "UPPMAX" or selection ==  "uppmax":
    selection = "bioinfo-tools CNVnator samtools vep bwa abyss fermikit vcftools"

print "creating conda environment"
FindSV_env="source activate FindSV_env\n"
command=["{} {} {}".format(os.path.join(programDirectory,"internal_scripts/CONDA/create_conda_env.sh"), os.path.join(programDirectory,"internal_scripts/CONDA/"), programDirectory)]
tmp=subprocess.check_output(command,shell = True)
    
if selection != "":
    FindSV_env += "module load {}\n".format(selection)
if thisroot != "":
    FindSV_env += "source {}".format(thisroot)

print "done! type source ./FindSV_env.sh before running FindSV to activate the FindSV environment, FindSV must be on during teh entire analysis, hence you might need to run a screen session to avoid interuption"

f= open("FindSV_env.sh", "w")
f.write(FindSV_env)
f.close()

print "installing nextflow"
os.system("curl -fsSL get.nextflow.io | bash")
os.system("chmod +x FindSV_env.sh")
