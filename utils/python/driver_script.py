# a python script to tie all of the CAPSEG steps together

# system level imports
import sys
import os
import argparse
import subprocess
import shutil

# create the command line arguments
parser = argparse.ArgumentParser()

# the analysis set name
parser.add_argument('-a', '--analysis_name',help="the analysis set name",required=True)

# bam files
parser.add_argument('-t', '--tumor',help="the tumor file list",required=True)
parser.add_argument('-n', '--normal',help="the normal file list",required=True)

# track information
parser.add_argument('-L', '--intervals',help="the intervals files",required=True)
parser.add_argument('-LB', '--intervals_bed',help="the intervals BED files",required=True)
parser.add_argument('-LC', '--intervals_csv',help="the intervals CSV file",required=True)

# reference information
parser.add_argument('-R', '--reference',help="the reference information",required=True)
parser.add_argument('-b', '--build',help="the genome build (hg18 or hg19)",required=True)

# other info  - directory info, etc
parser.add_argument('-o', '--output',help="the output file",required=True)
parser.add_argument('-temp', '--temp_dir',help="the temp directory",required=True)
parser.add_argument('-libdir', '--libdir',help="the library directory",required=True)
parser.add_argument('-queue', '--queue',help="the queue jar location",required=True)
parser.add_argument('-short_queue', '--short_queue',help="the temp directory",required=False,default="hour")
parser.add_argument('-uhd', '--use_historical_data',help="use the historical data (or not)",required=False,default="false")
parser.add_argument('-ktd', '--keep_temp_data',help="keep the temp data around",required=False,default="false")

args = parser.parse_args()

# some constants
job_project = "CAPSEG"
log_to_file = "logFile.txt"
HSL = "./"
iPerSplit = "400"
capseg_loc = args.libdir # + "/R/"
preprocessFile = "tumorData.txt"

# python location
py_location = "/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.6.5/bin/python"

# preprocess the data to transform the data into a table
preprocess_command = py_location + " " + args.libdir + "/utils/python/preprocess.py -t " + args.tumor + " -n " + args.normal + " -b " + args.build + " -o " + preprocessFile

print "About to run preprocessing with command: " + preprocess_command

p = subprocess.Popen(preprocess_command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

# make a local temp directory
local_tmp_dir = "./temp_directory"
try:
	os.mkdir(local_tmp_dir)
except:
	print "unable to create temp dir " + local_tmp_dir
	pass

# Check the return status
if sts != 0:
	print "Failed to process the data into a single table"

# are we generating coverage or are we going the distance?

# if we're good to go, start off the CAPSEG portion of the run
CAPSEG_command = "java -Xmx4g -Djava.io.tmpdir=" + args.temp_dir
CAPSEG_command += " -jar " + args.libdir + "/utils/Queue.jar -S "
CAPSEG_command += args.libdir + "/utils/CapSeg.scala -I " + preprocessFile
CAPSEG_command += " -L " + args.intervals + " -T " + args.intervals_bed
CAPSEG_command += " --tmpdir " + local_tmp_dir
CAPSEG_command += " -R " + args.reference
CAPSEG_command += " -bsub -build " + args.build + " -o ./ -jobQueue "
CAPSEG_command += args.short_queue + " -jobProject "
CAPSEG_command += job_project + " -memLimit 2 --log_to_file " + log_to_file
CAPSEG_command += " -iPerSplit 400 --longqueue " + args.queue
CAPSEG_command += " --baitcsv " + args.intervals_csv + " "
CAPSEG_command += " -retry 1 -asn "
CAPSEG_command += args.analysis_name + " -run"
CAPSEG_command += " -lib " + args.libdir + " -pl"
CAPSEG_command += " -uhd " + args.use_historical_data

print "About to run CAPSEG command: " + CAPSEG_command
#p = subprocess.Popen(CAPSEG_command, shell=True)
#sts = os.waitpid(p.pid, 0)[1]
p = subprocess.Popen(CAPSEG_command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

# remove the temp. dir we just setup
if not args.keep_temp_data.upper() == "TRUE":
        shutil.rmtree(local_tmp_dir)

local_q_dir = "./q_output"
try:
	os.mkdir(local_q_dir)
except:
	print "unable to create temp dir " + local_q_dir
	pass

# move all the Q script output to a seperate directory
for filename in os.listdir("."):
        if filename.startswith("Q") and (filename.endswith("txt") or filename.endswith("out")):
                shutil.move(filename, os.path.join(local_q_dir,filename))

print "returning CapSeg results of " + str(sts)
sys.exit(sts)
