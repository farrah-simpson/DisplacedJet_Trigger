#submit batch jobs for selection plots

import sys, os
import shutil
import getpass
import glob
import argparse
import subprocess
#set paths

parser = argparse.ArgumentParser(description='Submit batch jobs')
parser.add_argument('-i', "--input", default="", help="The input directory where the analyzer output trees are")
parser.add_argument('-o', "--output", default="", help="The output directory for selection plots")
parser.add_argument('-d', "--data", default=False, action="store_true", help="data/MC to determine the name of triggers")
parser.add_argument("-a","--abs", default=False, action="store_true", help="use absolute path")
parser.add_argument("-w", "--weight", default="1.0",  help="The weighte factor of the sample")
 
args=parser.parse_args()
weight = args.weight


current=os.getcwd()
basefolder=args.input
bashjob="trigeff.csh"
pathbashjob="{0}/{1}".format(current, bashjob)
pyscript="Trigger_trk_AN.py"
pathpyscript="{0}/{1}".format(current, pyscript)

root_files = glob.glob("{0}/*.root".format(basefolder))
if args.abs:
    root_files=glob.glob("/eos/uscms{0}/*.root".format(basefolder))
    #root_files = subprocess.check_output(["eosls", "{0}/*.root".format(basefolder)])

for rootfile in root_files:
    print("+look at the root file: {0}".format(rootfile))
    if args.abs:
        rootfile=rootfile[rootfile.find("/store"): ]
        rootfile="root://cmseos.fnal.gov/"+rootfile
    #rootfile = rootfile.lstrip("/eos/uscms")
    folder0 = rootfile.split("/")[-1].rstrip(".root")
    folder1 = rootfile.split("/")[-2]

    folder = args.output+"/"
    if not os.path.isdir(folder): os.mkdir(folder)
    folder = folder+folder0
    if not os.path.isdir(folder): os.mkdir(folder)

    os.chdir(folder)

    shutil.copyfile(pathpyscript, pyscript)
    shutil.copyfile(pathbashjob, bashjob)

    if args.data:
        DATA="-d"
    else:
        DATA=""
    condor_filename="analyze_condor_{0}".format(folder0)
    fcondor=open(condor_filename, "w")
    fcondor.write("Executable = {0}\n".format(bashjob))
    fcondor.write("Universe = vanilla\n")
    fcondor.write("transfer_input_files = {0}\n".format(pyscript)) 
    fcondor.write("should_transfer_files = YES\n")
    fcondor.write("Output = {0}/{1}/run.out\n".format(current, folder))
    fcondor.write("Error  = {0}/{1}/run.err\n".format(current, folder))
    fcondor.write("Log    = {0}/{1}/run.log\n".format(current, folder))
    if args.abs:
        #rootfile="root://cmsxrootd.fnal.gov/"+rootfile
        fcondor.write("Arguments = {0} {1} {2} {3} {4}\n".format(pyscript, rootfile, folder1+"_"+folder0,weight, DATA))
    else:
        fcondor.write("Arguments = {0} {1} {2} {3} {4}\n".format(pyscript, current+"/"+rootfile, folder1, weight, DATA))
    fcondor.write("Queue\n")
    fcondor.close()
    
    os.system("chmod +x {0} {1} analyze_condor_{2}".format(bashjob, pyscript, folder0))
    os.system("condor_submit analyze_condor_{0}".format(folder0))
 
    os.chdir(current)
