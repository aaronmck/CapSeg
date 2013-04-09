# ------------------------------------------------------------------------
# take the input VCF, and split it out (using the specified sample list
# mapping) into individual BED files
# ------------------------------------------------------------------------

# system level imports
import sys
import os
import argparse
import shutil

# create the command line arguments
parser = argparse.ArgumentParser()

# the mapping of sample name to BAM sample name; they are not always identical, and we want a way to map them into CapSeg expected names
parser.add_argument('-map', '--mapping', help="the mapping of samples we want to their BAM sample name (not always the same)", required=True)

# where to dump the output: we use the sample name plus .bed
parser.add_argument('-out', '--outputdir', help="where to put the output files", required=True)

# our VCF file
parser.add_argument('-vcf', '--vcf', help="the input vcf file", required=True)

args = parser.parse_args()

# find the header line in the input VCF file
vcf_input = open(args.vcf)

# load up the provided mapping file - map the BAM sample name to the BAM file
sample_mapping = {}
sm = open(args.mapping,"r")
hdr = sm.readline()
for line in sm:
    sp = line.strip().split("\t")
    sample_mapping[sp[1]] = sp[0]

sample_to_position = {}
sample_to_output = {}

header_loaded = False
for line in vcf_input:
    if line.startswith("##"):
        pass
    elif line.startswith("#"):
        sp = line.strip().split("\t")
        # looking for the following line #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1....
        if len(sp) < 10:
            raise NameError("Unable to process the VCF: it has only " + str(len(sp)) + " columns, we expect at least 10")

        # now find each of the samples
        samples = sp[9:len(sp)]
        for key, value in sample_mapping.iteritems():
            if not key in samples:
                raise NameError("Sample name " + key + " missing from the input file")
            sample_to_position[key] = samples.index(key)
        header_loaded = True
        
        for key, value in sample_mapping.iteritems():
            flname = open(os.path.join(args.outputdir,key+".bed"),"w")
            sample_to_output[key] = flname

    # process the rest of the lines in the file        
    else:
        if header_loaded != True:
            raise NameError("Header not seen in the file")
        # check that the site is a PASS
        sp = line.strip().split("\t")
        if sp[6] != "PASS":
            continue

        # get the position
        chrom = sp[0]
        pos = sp[1]
        
        # now check each sample; if they're a het print it to the respective file
        for key, val in sample_to_position.iteritems():
            geno_full = sp[sample_to_position[key] + 9]
            geno_each = geno_full.split(":")[0].split("/")
            if len(geno_each) != 2 or geno_each[0] == geno_each[1]:
                continue

            # ok, we have a het to output, dump it to the target file
            output_file = sample_to_output[key]
            output_file.write(chrom + "\t" + str(int(pos)-1) + "\t" + pos + "\t" + (geno_full) + "\n")
            
        

