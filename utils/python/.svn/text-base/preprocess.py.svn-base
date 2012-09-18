# ------------------------------------------------------------------------------------------
# preprocess the tumor and normal bam files, and try to extract the samples from each file
# matching up tumors and normals
# ------------------------------------------------------------------------------------------

# system level imports
import sys
import os
import argparse
import subprocess

# local imports

# create the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tumor',help="the tumor file list",required=True)
parser.add_argument('-n', '--normal',help="the normal file list",required=True)
parser.add_argument('-alt', '--alternate_normals',help="the normal file report",required=False)
parser.add_argument('-b', '--build',help="the genome build (hg18 or hg19)",required=False)
parser.add_argument('-o', '--output',help="the output file",required=False)

args = parser.parse_args()

# run the command:
# samtools view -H /seq/picard_aggregation/C335/CW84S/v2/CW84S.bam | grep "@RG"
normals = open(args.normal,"r")
tumors = open(args.tumor,"r")

# read in the tumor and normal header line
#normal_header = normals.readline()
#tumor_header = tumors.readline()

nmap = {}
tmap = {}

def getSamples(file):
	samples = {}
	p = os.popen('samtools view -H ' + file + ' | grep "@RG"',"r")
	while 1:
		line = p.readline()
		if not line: break
		sp = line.split()
		for token in sp:
			if token.startswith("SM"):
				samples[token.split(":")[1]] = file
	return samples
			
			
normal_bam_files = []
tumor_bam_files = []

for line in normals:
	sp = line.strip().split()
	if len(sp) > 1 and os.path.exists(sp[1]):
		sm = sp[0].rstrip("-Normal")
		if nmap.has_key(sm) or sp[1] in normal_bam_files:
			print "repeated normal " + sm + " bam " + sp[1]
		else:
			nmap[sm] = sp[1]
			normal_bam_files.append(sp[1])
	elif len(sp) > 1:
		print "Doesn't exist " + sp[1]
	
for line in tumors:
	sp = line.strip().split()
	if len(sp) > 1 and os.path.exists(sp[1]):
		sm = sp[0].rstrip("-Tumor")
		if tmap.has_key(sm) or sp[1] in tumor_bam_files:
			print "repeated tumor " + sm + " bam " + sp[1]
		else:
			tmap[sm] = sp[1]
			tumor_bam_files.append(sp[1])
	elif len(sp) > 1:
		print "Doesn't exist " + sp[1]
		
samples = set()
for key,value in tmap.iteritems():
	samples.add(key)
for key,value in nmap.iteritems():
	samples.add(key)

output = open(args.output,"w")
output.write("sample\ttumor_bam\tnormal_bam\n")
for sample in samples:
	if nmap.has_key(sample) and tmap.has_key(sample):
		output.write(sample + "\t" + tmap[sample] + "\t" + nmap[sample] + "\n")
	elif tmap.has_key(sample):
		output.write(sample + "\t" + tmap[sample] + "\tNA\n")
	else: 
		output.write(sample + "\tNA\t" + nmap[sample] + "\n")
output.close()
