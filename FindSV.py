import argparse
import sys
import yaml
import os
import re
import fnmatch
import subprocess

def get_output_dir(config):
    working_dir = None
    for line in open(config):
        param=line.split("=")
        if "working_dir" in param[0]:
            working_dir= re.split('\"|\'',param[1])[1]
        
    return working_dir
    
def retrieve_status(status,nextflow_output):
    for sample in nextflow_output:
        sample_id=sample.split("\n")[0].replace(".bam","").split("/")[-1]
        if "Error executing process" in sample:
            status[sample_id]["status"]="FAILED:CALLING"
            if "Submitted process > annotate" in sample:
                 status[sample_id]["status"]="FAILED:ANNOTATION"
        else:
            status[sample_id]["status"]="SUCCESS"
        status[sample_id]["message"]=sample
        
    return status
    
def worker(bam_files,args,status):
    print "Processing, please do not turn FindSV off"
    
    processes=[]
    for sample in bam_files:
        prefix=sample["bam"].split("/")[-1].replace(".bam","")
        status[prefix]={}
        status[prefix]["bam_file"]=sample["bam"]
        f = os.tmpfile()
        e = os.tmpfile()
        if sample["mode"] == "full":
            process=["./launch_core.sh",sample["bam"],args.config,args.output]
        elif sample["mode"] == "annotate":
            process=["./launch_core_reannotate.sh",sample["bam"],args.config,args.output,os.path.join(args.output,prefix+"_CombinedCalls.vcf")]
        p = subprocess.Popen(process,stdout=f,stderr=e)
        processes.append((p, f,e))
            
    nextflow_output=[]
    for p, f ,e in processes:
        p.wait()
        f.seek(0)
        nextflow_output.append(f.read())
        f.close()
    status.update( retrieve_status(status,nextflow_output) )
    
    return status

def print_yaml(status,working_dir):
    for sample in status:
        status[sample]["combined_caller_vcf"]=os.path.join(working_dir,sample+"_CombinedCalls.vcf")
        status[sample]["annotated_vcf"]=os.path.join(working_dir,sample+"_FindSV.vcf")
        if "FAILED" in status[sample]["status"]:
            print sample
            print status[sample]["status"]
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    f = open(os.path.join(working_dir,"tracker.yml"), 'w')
    f.write(yaml.dump(status))
    f.close()
    print "DONE"

def read_yaml(working_dir):
    status={}
    if os.path.exists( os.path.join(working_dir,"tracker.yml") ):
        with open(os.path.join(working_dir,"tracker.yml"), 'r') as stream:
            status=yaml.load(stream)
    return status

parser = argparse.ArgumentParser("FindSV core module",add_help=False)
parser.add_argument('--bam', type=str,help="analyse the bam file using FindSV")
parser.add_argument("--folder", type=str,help="analyse every bam file within a folder using FindSV")
parser.add_argument('--output', type=str,default=None,help="the output is stored in this folder, default output folder is fetched from the config file") 
parser.add_argument("--config",type=str,required=True, default=None,help="the config file")
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
    
    status=read_yaml(args.output)
    prefix=args.bam.split("/")[-1].replace(".bam","")
    if not prefix in status:
        status=worker([{"bam":args.bam,"mode":"full"}],args,status)
        print_yaml(status,args.output)
    else:
        print "Error: the bam file is already analysed, either restart the analysis using restart, or use an other working directory"
        
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
    
    processes = []
    bam_files=[]
    for root, dirnames, filenames in os.walk(args.folder):
        for filename in fnmatch.filter(filenames, '*.bam'):
            bam_file=os.path.join(root, filename)
            prefix=bam_file.split("/")[-1].replace(".bam","")
            if not prefix in status:
                bam_files.append({"bam":bam_file,"mode":"full"})
    
    if bam_files:      
        status=worker(bam_files,args,status)
        print_yaml(status,args.output)
    else:
        print "error: no new bam files was found, use the restart module if you wish to restart any sample"
        
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
    
    #restart everything
    if args.full:
        bam_files=[]
        for sample in status:
            bam_files.append({"bam":status[sample]["bam_file"],"mode":"full"})
        
        status={}
        if bam_files:      
            status=worker(bam_files,args,status)
            print_yaml(status,args.output)
            
    elif args.annotation:
        bam_files=[]
        for sample in status:
            bam_files.append({"bam":status[sample]["bam_file"],"mode":"annotate"})
        
        status={}
        if bam_files:      
            status=worker(bam_files,args,status)
            print_yaml(status,args.output)    
                     
    elif args.failed:
    
        bam_files=[]
        for sample in status:
            if "FAILED:ANNOTATION" in status[sample]["status"]:
               bam_files.append({"bam":status[sample]["bam_file"],"mode":"annotate"})
            elif "FAILED:CALLING" in status[sample]["status"]:
                bam_files.append({"bam":status[sample]["bam_file"],"mode":"full"})
        
        if bam_files:      
            status=worker(bam_files,args,status)
            print_yaml(status,args.output) 
            
    else:
        parser.print_help()
                
else:
    parser.print_help()
