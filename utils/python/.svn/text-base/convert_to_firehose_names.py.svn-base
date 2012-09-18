# a python script to create load files from the CapSeg output data

# system level imports
import sys
import os
import argparse
import subprocess
import shutil
# create the command line arguments
parser = argparse.ArgumentParser()

# the analysis set name
parser.add_argument('-m', '--master_file',help="the firehose sample name to tumor and normal bam file mapping",required=True)
parser.add_argument('-b', '--bam_to_sample',help="the bam file to sample name mapping",required=True)
parser.add_argument('-s', '--seg_file_location',help="where to find the segmentation results",required=True)
parser.add_argument('-o', '--output',help="the output file",required=True)

args = parser.parse_args()

# open the master file, and get the mapping of individuals to tumor bam files
master_fl = open(args.master_file,"r")
header = master_fl.readline()
if header.strip() != "sample\ttumor_bam\tnormal_bam":
    raise NameError("Unable to match appropriate header in master file, we expect: sample<tab>tumor_bam<tab>normal_bam")

tumor_to_ind = {}
for line in master_fl:
    sp = line.strip().split("\t")
    tumor_to_ind[os.path.basename(sp[1])] = sp[0]

# now get the mapping of tumor bam file to sample name
bam_to_sample_fl = open(args.bam_to_sample,"r")
header = bam_to_sample_fl.readline()
if header.strip() != "BAM\tSample":
    raise NameError("Unable to match appropriate header in bam to sample file, we expect: BAM<tab>Sample")

sample_to_bam = {}
for line in bam_to_sample_fl:
    sp = line.strip().split("\t")
    sp[1] = sp[1].replace(" ",".")
    sp[1] = sp[1].replace("-",".")
    sample_to_bam[sp[1]] = sp[0]

# now make a mapping of an individual to it's tumors segmentation file
final_mapping = {}
for fl in os.listdir(args.seg_file_location):
    tm_name = fl.strip()
    tm_name = tm_name[0:len(tm_name)-8]
    if not sample_to_bam.has_key(tm_name):
        tm_name = tm_name[1:len(tm_name)]
    if not sample_to_bam.has_key(tm_name):
        print "dropping seg file " + fl + " (" + tm_name + ") since it's not in the sample_to_bam mapping"
        continue
    bam = sample_to_bam[tm_name]
    if not tumor_to_ind.has_key(bam):
        print "dropping seg file " + fl + " since it's not in the bam to individual mapping"
        continue
    ind = tumor_to_ind[bam]
    final_mapping[ind] = os.path.join(args.seg_file_location,fl)


# now make the output file to load into filehose
output = open(args.output,"w")
output.write("individual_id\tcapseg_segmentation_file\n")
for key,val in final_mapping.iteritems():
    output.write(key + "\t" + os.path.abspath(val) + "\n")


