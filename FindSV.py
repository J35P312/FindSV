import argparse
import sys
import yaml
import os
import re
import fnmatch
import subprocess
import tempfile
import glob

def get_output_dir(config):
    working_dir = None
    for line in open(config):
        param=line.split("=")
        if "working_dir" in param[0]:
            working_dir= re.split('\"|\'',param[1])[1]
        
    return working_dir

def retrieve_bam(working_dir):
    programDirectory = os.path.dirname(os.path.abspath(__file__))
    TIDDIT=os.path.join(programDirectory,"work",working_dir)+"**/*.bam"
    bam=glob.glob(os.path.join(TIDDIT))[0]
    prefix=bam.split("/")[-1].split(".bam")[0]
    return([prefix,bam])
    
def retrieve_status(status,nextflow_output):

    tmp_dict={}
    process_2_prefix={}
    for message in nextflow_output:
        for line in message.split("\n"):
            if "Submitted" in line:
                content=line.strip().split()
                prefix,bam=retrieve_bam(content[0][1:-1])
                if not prefix in tmp_dict:
                    tmp_dict[prefix]={}
                    tmp_dict[prefix]["status"]="SUCCESS"
                tmp_dict[prefix][content[-2]]={}
                tmp_dict[prefix][content[-2]]["working_dir"]=content[0][1:-1]
                tmp_dict[prefix][content[-2]]["status"]="SUCCESS"
                tmp_dict[prefix][content[-2]]["message"]=line.strip()
                
                if not content[-2] in process_2_prefix:
                    process_2_prefix[content[-2]] = {}
                process_2_prefix[content[-2]][content[-1]]=prefix
                    
                
            elif "terminated with an error" in line:
                content=line.split("`")[1].split()   
                
                process= "CALLING"
                if content[1] == "annotate":
                    process="ANNOTATION"  
                prefix=process_2_prefix[content[0]][content[1]]
                
                tmp_dict[prefix][content[0]]["status"]="FAILED:{}".format(process)
                tmp_dict[prefix][content[0]]["message"] += "\n" + line.strip()
                tmp_dict[prefix]["status"]="FAILED:{}".format(process)
            else:
                print "FAILURE"
                print message.strip()
                break
            
    for sample in tmp_dict:
            	        	
        sample_id,bam=retrieve_bam(tmp_dict[sample]["TIDDIT"]["working_dir"])
        #we are reruning samples, update the sample info
        if sample_id in status:
            for step in tmp_dict[sample]:
                    status[sample_id][step]=tmp_dict[sample][step]
        #we are analysing new sample, create the sample info
        else:
            status[sample_id]=tmp_dict[sample]
    return status
    
def worker(bam_files,args,status):
    print "Processing, please do not turn FindSV off"
    nextflow_path = os.path.join( os.path.dirname(os.path.abspath(__file__)),"nextflow")
    
    processes=[]
    for sample in bam_files:
        f = tempfile.NamedTemporaryFile()
        e = tempfile.NamedTemporaryFile()
        if sample["mode"] == "full":
            process=["./launch_core.sh",sample["bam"],args.config,args.output,nextflow_path]
        elif sample["mode"] == "annotate":
            process=["./launch_core_reannotate.sh",sample["bam"],args.config,args.output,sample["vcf"],nextflow_path]
        p = subprocess.Popen(process,stdout=f,stderr=e)
        processes.append((p, f,e))
            
    nextflow_output=[]
    for p, f ,e in processes:
        try:
            print ("waiting for {} to finish... please hold".format(f.name))
            p.wait()
        except:
            print "FAIL, unnable to contact:{}".format(f.name)
        f.seek(0)
        nextflow_output.append(f.read())
        f.close()
    status =  retrieve_status(status,nextflow_output)
    
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
        status[prefix]={}
        status[prefix]["bam_file"]=args.bam
        status[prefix]["status"]="SUBMITTED"
        print_yaml(status,args.output)
    
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
    
        print_yaml(status,args.output)
        bam_files=",".join(bam_files) 
        status=worker([{"bam":bam_files,"mode":"full"}],args,status)
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
            bam_files.append(status[sample]["bam_file"])

       
        status={}
        for bam in bam_files:
            prefix=bam.split("/")[-1].replace(".bam","")
            status[prefix]={}
            status[prefix]["bam_file"]=bam
            status[prefix]["status"]="SUBMITTED" 
     
        if bam_files: 
            print_yaml(status,args.output)
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
            print_yaml(status,args.output)
            bam_files=",".join(bam_files)
            vcf_files=",".join(vcf_files)
            status=worker([{"bam":bam_files,"vcf":vcf_files,"mode":"annotate"}],args,status)
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
            elif "FAILED" in status[sample]["status"]:
               full.append(status[sample]["bam_file"])
               status[sample]["status"]="SUBMITTED"
        
        bam_files=[]
        if full:
            bam_files.append({"bam":",".join(full),"mode":"full"})
        if annotation_bam:
            bam_files.append({"bam":",".join(fannotation_bam),"vcf":",".join(fannotation_vcf),"mode":"annotate"})
        if bam_files:
            print_yaml(status,args.output)    
            status=worker(bam_files,args,status)
            print_yaml(status,args.output) 
            
    else:
        parser.print_help()
                
else:
    parser.print_help()
