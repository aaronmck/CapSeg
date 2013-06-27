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
parser.add_argument('-v', '--vcf',help="the vcf file to pull a sample listing from",required=True)
parser.add_argument('-alt', '--alternate_normals',help="the normal file report",required=False)
parser.add_argument('-b', '--build',help="the genome build (hg18 or hg19)",required=False)
parser.add_argument('-o', '--output',help="the output file",required=False)

args = parser.parse_args()

# run the command:
# samtools view -H /seq/picard_aggregation/C335/CW84S/v2/CW84S.bam | grep "@RG"
normals = open(args.normal,"r")
tumors = open(args.tumor,"r")

# read in the tumor and normal header line
nmap = {}
tmap = {}


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# helper function section
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
def getBAMSamples(input_file):
        """Gets sample names from the header of the bam file provided


        Args:
            input_file: the BAM input file to load up the sample names from

        Returns:
            a list of sample names found, or if no samples are found or the file doesn't exist, None
        """
	samples = {}
        if not os.path.exists(input_file):
                return None
	p = os.popen('samtools view -H ' + input_file + ' | grep "@RG"',"r")
	while 1:
		line = p.readline()
		if not line: break
		sp = line.split()
		for token in sp:
			if token.startswith("SM"):
				samples[token.split(":")[1]] = file
	return samples

def getVCFSamples(input_file):
        """Gets sample names from the VCF file they provide

        Args:
            input_file: the VCF input file to load up the sample names from

        Returns:
            a list of sample names found, or if no samples are found or the file doesn't exist, None
        """
	samples = {}
        if not os.path.exists(input_file):
                return None
	p = open(input_file,"r")
        for line in p:
		if not line.startswith("#"):
                       return None
                if line.startswith("#CHROM"):
                        sp = line.strip().split("\t")
                        if len(sp) < 10:
                                return None
                        return sp[9:len(sp)]
        return None


normal_bam_files = []
tumor_bam_files = []

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# load up the tumors and the normals from the input lists from Firehose
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#
# now go and get the samples from each bam file, and check if it's in the VCF file or not
# this is an annoying check, but worth it to determine that we're only processing samples
# that have calls.
#
# also collect a list of all the samples we've seen in either the tumor or the normal
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
normal_sample_to_vcf = {}
vcf_samples = getVCFSamples(args.vcf)

samples = set()
for key,value in tmap.iteritems():
	samples.add(key)
for key,value in nmap.iteritems():
	samples.add(key)
        # get the normal sample from the bam file, if it exists
        normal_samples = getBAMSamples(value)
        if len(normal_samples) != 1:
                print("unable to find an appropriate number of samples from " + value + "; dropping")
        else:
                # check that the sample is in the VCF
                if normal_samples[0] in vcf_samples:
                        normal_sample_to_vcf[key] = args.vcf





output = open(args.output,"w")
output.write("sample\ttumor_bam\tnormal_bam\tvcf_file\n")
for sample in samples:
        output_str = sample + "\t"

        # add the tumor if it exists
	if tmap.has_key(sample):
		output_str += tmap[sample] + "\t"
        else:
                output_str += "NA\t"

        # add the normal if it exists
	if nmap.has_key(sample):
		output_str += nmap[sample] + "\t"
        else:
                output_str += "NA\t"

        # add the vcf if it exists
	if normal_sample_to_vcf.has_key(sample):
		output_str += normal_sample_to_vcf[sample] + "\n"
        else:
                output_str += "NA\n"
        output.write(output_str)

output.close()
