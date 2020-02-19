import argparse
import sys
import yaml
import os
import re
import fnmatch
import glob

def get_output_dir(config):
    working_dir = None
    for line in open(config):
        param=line.split("=")
        if "working_dir" in param[0]:
            working_dir= re.split('\"|\'',param[1])[1]
        
    return working_dir
    
def retrieve_status(status,trace,working_dir):

    tmp_dict={}
    process_2_prefix={}
    first = True
    for line in open(trace):
        if first:
            first = False
            continue  
        content=line.split("\t")
        bam=content[-1]
        prefix=bam.split("/")[-1].split(".bam")[0]
        task=content[2].split()[0]
        
        if not prefix in tmp_dict:
            tmp_dict[prefix]={}
            tmp_dict[prefix]["status"]="SUCCESS"
        tmp_dict[prefix][task]={}
        tmp_dict[prefix][task]["working_dir"]=os.path.join(content[1])
        tmp_dict[prefix][task]["status"]=content[-2]
            
    for sample in tmp_dict:
        if "Assemblatron" in tmp_dict[sample]:
            if  not tmp_dict[sample]["Assemblatron"]["status"] == "COMPLETED":
                tmp_dict[sample]["status"] = "FAILED:CALLING"
			
        if "TIDDIT" in tmp_dict[sample]:
            if  not tmp_dict[sample]["TIDDIT"]["status"] == "COMPLETED":
                tmp_dict[sample]["status"] = "FAILED:CALLING"
                
        if "CNVnator" in tmp_dict[sample]:
            if  not tmp_dict[sample]["CNVnator"]["status"] == "COMPLETED":
                tmp_dict[sample]["status"] = "FAILED:CALLING"
                
        if "combine" in tmp_dict[sample]:
            if  not tmp_dict[sample]["combine"]["status"] == "COMPLETED":
                tmp_dict[sample]["status"] = "FAILED:CALLING"
        if "annotate" in tmp_dict[sample]:
            if  not tmp_dict[sample]["annotate"]["status"] == "COMPLETED":
                tmp_dict[sample]["status"] = "FAILED:ANNOTATION"  

        if sample in status:
            for step in tmp_dict[sample]:
                    status[sample][step]=tmp_dict[sample][step]
        #we are analysing new sample, create the sample info
        else:
            status[sample]=tmp_dict[sample]
                     
    return status
    
def worker(bam_files,args,status):
    path=os.path.dirname(sys.argv[0])
    if not path:
        path = "."
    print ("Processing, please do not turn off FindSV")

    if os.path.isfile("{}/trace.txt".format(args.output)):
        os.remove("{}/trace.txt".format(args.output))
    if os.path.isfile("{}/log_annotation.txt".format(args.output)):
        os.remove("{}/log_annotation.txt".format(args.output))
    if os.path.isfile("{}/trace_call.txt".format(args.output)):
        os.remove("{}/trace_call.txt".format(args.output))

    os.system("mkdir {}".format(args.output))
    f=open("{}/lock".format(args.output),"w")
    f.write("If this file is not deleted, the analysis is either ongoing, or it crashed.\n You need to remove this file in order to retry the analysis.\n Do not delete this file if the analyis is ongoing!")
    f.close()
    
    sample=bam_files[0]
    if sample["mode"] == "full":
        process=["{}/launch_core.sh".format(path),sample["bam"],args.config,args.output]
    elif sample["mode"] == "annotate":
        process=["{}/launch_core_reannotate.sh".format(path),sample["bam"],args.config,args.output,sample["vcf"]]
    elif sample["mode"] == "restart_failed":
        process=["{}/launch_core_restart_failed.sh".format(path),sample["bam_call"],sample["bam_annotate"],args.config,args.output,sample["vcf"]]
    os.system(" ".join(process))
    
    status =  retrieve_status(status, os.path.join(args.output,"trace.txt"),args.output)
    return status

def print_yaml(status,working_dir):
    for sample in status:
        status[sample]["combined_caller_vcf"]=os.path.join(working_dir,sample+"_CombinedCalls.vcf")
        status[sample]["annotated_vcf"]=os.path.join(working_dir,sample+"_FindSV.vcf")
        if "FAILED" in status[sample]["status"]:
            print (sample)
            print (status[sample]["status"])

    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    if os.path.isfile(os.path.join(working_dir,"tracker.yml")):
        os.system("cp {} {}".format(os.path.join(working_dir,"tracker.yml"),os.path.join(working_dir,"tracker.bkp")))
    f = open(os.path.join(working_dir,"tracker.yml"), 'w')
    f.write(yaml.dump(status))
    f.close()
    os.remove("{}/lock".format(working_dir))

    print ("DONE")

def check_lock(working_dir):
    if os.path.isfile(os.path.join(working_dir,"lock")):
        print ("error: the ouput directory is locked! Choose another output directory, or wait for the ongoing analysis to finish")
        print ("If the analysis crashed, remove the lock file and restart the analysis")
        quit()

def read_yaml(working_dir):
    status={}
    if os.path.exists( os.path.join(working_dir,"tracker.yml") ):
        with open(os.path.join(working_dir,"tracker.yml"), 'r') as stream:
            status=yaml.load(stream)
    return status

parser = argparse.ArgumentParser("FindSV",add_help=False)
parser.add_argument('--bam', type=str,help="analyse the bam file using FindSV")
parser.add_argument("--folder", type=str,help="analyse every bam file within a folder using FindSV")
parser.add_argument('--output', type=str,default=None,help="the output is stored in this folder, default output folder is fetched from the config file") 
parser.add_argument("--config",type=str, default=None,help="the config file")
parser.add_argument("--restart",action="store_true",help="restart module: perform the selected restart on the specified folder")
args, unknown = parser.parse_known_args()
programDirectory = os.path.dirname(os.path.abspath(__file__)) 
#analyse one single bam file



if args.bam:
    parser = argparse.ArgumentParser("FindSV core module",add_help=False)
    parser.add_argument('--bam', type=str,help="analyse the bam file using FindSV")
    parser.add_argument('--output', type=str,default=None,help="the output is stored in this folder")
    parser.add_argument("--config",type=str, required=True,default=None,help="the location of the config file(default= the same folder as the FindSV-core script")
    args = parser.parse_args()
    
    if not args.output:
        args.output=get_output_dir(args.config)
    check_lock(args.output)
    
    status=read_yaml(args.output)
    prefix=args.bam.split("/")[-1].replace(".bam","")
    if not prefix in status:
        status[prefix]={}
        status[prefix]["bam_file"]=args.bam
        status[prefix]["status"]="SUBMITTED"
    
        status=worker([{"bam":args.bam,"mode":"full"}],args,status)
        print_yaml(status,args.output)
    else:
        print ("Error: the bam file is already analysed, either restart the analysis using restart, or use an other working directory")
        
#analyse all bamfiles within a folder(recursive searching)
elif args.folder:
    parser = argparse.ArgumentParser("FindSV core module",add_help=False)
    parser.add_argument("--folder", type=str,help="analyse every bam file within a folder using FindSV")
    parser.add_argument('--output', type=str,default=None,help="the output is stored in this folder")
    parser.add_argument("--config", required=True, type=str, default=None,help="the config file")
    args = parser.parse_args()
    
    if not args.output:
        args.output=get_output_dir(args.config)
    status=read_yaml(args.output)
    check_lock(args.output)

    bam_files=[]
    for root, dirnames, filenames in os.walk(args.folder):
        for filename in fnmatch.filter(filenames, '*.bam'):
            bam=os.path.join(root, filename)
            prefix=bam.split("/")[-1].replace(".bam","")
            if not prefix in status:
                status[prefix]={}
                status[prefix]["bam_file"]=bam
                status[prefix]["status"]="SUBMITTED"
                bam_files.append(bam)
    
    if bam_files: 
    
        bam_files=",".join(bam_files) 
        status=worker([{"bam":bam_files,"mode":"full"}],args,status)
        print_yaml(status,args.output)
        
    else:
        print ("error: no new bam files was found, use the restart module if you wish to restart any sample")
        
#the restart module
elif args.restart:
    parser = argparse.ArgumentParser("FindSV core module:restart module")
    parser.add_argument("--failed",action="store_true",help="restart all failed samples")
    parser.add_argument("--annotation",action="store_true",help="reruns the annotation step on all samples")
    parser.add_argument("--full",action="store_true",help="restarts the analysis from scratch")
    parser.add_argument("--restart",type=str,nargs="*",help="restart module: perform the selected restart on the specified folder(default output if none is chosen)")
    parser.add_argument('--output', type=str,default=None,help="the output is stored in this folder")
    parser.add_argument("--config",type=str, default=None,help="the config file")
    args = parser.parse_args()

    if not args.output:
        args.output=get_output_dir(args.config)
    status=read_yaml(args.output)
    check_lock(args.output)
    
    #restart everything
    if args.full:
        bam_files=[]
        for sample in status:
            bam_files.append(status[sample]["bam_file"])

       
        status={}
        for bam in bam_files:
            prefix=bam.split("/")[-1].replace(".bam","")
            status[prefix]={}
            status[prefix]["bam_file"]=bam
            status[prefix]["status"]="SUBMITTED" 
     
        if bam_files: 
            bam_files=",".join(bam_files)     
            status=worker([{"bam":bam_files,"mode":"full"}],args,status)
            print_yaml(status,args.output)
            
    elif args.annotation:
        bam_files=[]
        vcf_files=[]
        for sample in status:
            if status[sample]["status"] == "SUCCESS" or status[sample]["status"] == "FAILED:ANNOTATION":
                bam_files.append(status[sample]["bam_file"])
                vcf_files.append(status[sample]["combined_caller_vcf"])
                status[sample]["status"]="SUBMITTED"
                
        if bam_files:
            bam_files=",".join(bam_files)
            vcf_files=",".join(vcf_files)
            status=worker([{"bam":bam_files,"vcf":args.output,"mode":"annotate"}],args,status)
            print_yaml(status,args.output)    
                     
    elif args.failed:
    
        full=[]
        annotation_bam=[]
        annotation_vcf=[]
        for sample in status:
            if "FAILED:ANNOTATION" in status[sample]["status"]:
                annotation_bam.append(status[sample]["bam_file"])
                annotation_vcf.append(status[sample]["combined_caller_vcf"])
                status[sample]["status"]="SUBMITTED"
            if "FAILED:CALLING" in status[sample]["status"]:
               full.append(status[sample]["bam_file"])
               status[sample]["status"]="SUBMITTED"
        
        bam_files=[]
        if full or annotation_bam:
            bam_files= [ {"bam_call":",".join(full),"bam_annotate":",".join(annotation_bam),"vcf":args.output,"mode":"restart_failed"} ]
            if not full:
                bam_files[0]["bam_call"] = "NONE"
            if not annotation_bam:
                bam_files[0]["bam_annotate"] = "NONE"
                bam_files[0]["vcf"] = "NONE"
  
            status=worker(bam_files,args,status)
            print_yaml(status,args.output) 
            
    else:
        parser.print_help()
                
else:
    parser.print_help()
