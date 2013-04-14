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

# the mapping of normal sample name to BAM sample name; they are not always identical, and we want a way to map them into CapSeg expected names
parser.add_argument('-nmap', '--normal_mapping', help="the mapping of normal samples to their BAM sample name (not always the same)", required=True)

# the mapping of the bam files of tumors and normals to their individual name; we need this to tie the normals and tumors together
parser.add_argument('-imap', '--individual_mapping', help="individual mapping to their tumor and normal bam files", required=True)

# the mapping of the tumor sample name to BAM sample name; they are not always identical, and we want a way to map them into CapSeg expected names
parser.add_argument('-tmap', '--tumor_mapping', help="the mapping of tumor samples to their BAM sample name (not always the same)", required=True)

# where to dump the output: we use the sample name plus .bed
parser.add_argument('-out', '--outputdir', help="where to put the output files", required=True)

# our VCF file
parser.add_argument('-vcf', '--vcf', help="the input vcf file", required=True)

args = parser.parse_args()

# find the header line in the input VCF file
vcf_input = open(args.vcf)

# --------------------------------------------------------------------------------------------------------------------------------------------------------
# what we'd like to do is map the normal sample names, what we'll see in the VCF/BAM, to their tumor sample names, also from the BAM/VCF file. This means going from
# the bam -> sample mapping, finding the bam, finding the bam in the ind -> bam list, and going back in the tumor mapping and finding that sample
# --------------------------------------------------------------------------------------------------------------------------------------------------------

# the actual final mapping -- we'll layer on each relation in turn
normal_sample_to_bam = {}
normal_bam_to_tumor_bam = {}
tumor_bam_to_sample = {}

final_normal_sample_to_tumor_sample = {}

# simple case, sample to bam in normal
sm = open(args.normal_mapping,"r")
hdr = sm.readline()
for line in sm:
    sp = line.strip().split("\t")
    normal_sample_to_bam[sp[1]] = sp[0]

# now go from normal bam to tumor bam
sm = open(args.individual_mapping,"r")
hdr = sm.readline()
for line in sm:
    sp = line.strip().split("\t")
    normal_bam_to_tumor_bam[sp[2]] = sp[1]

# final go from tumor bam to tumor sample
sm = open(args.tumor_mapping,"r")
hdr = sm.readline()
for line in sm:
    sp = line.strip().split("\t")
    tumor_bam_to_sample[sp[0]] = sp[1]

# now create the final mapping for each normal sample we have
for key,value in normal_sample_to_bam:
    if not normal_bam_to_tumor_bam.has_key(value):
        raise NameError("Unable to map the normal sample " + key + " to a tumor bam file")
    tumor_bam = normal_bam_to_tumor_bam[key]
    
    if not tumor_bam_to_sample.has_key(value):
        raise NameError("Unable to map the tumor bam " + tumor_bam + " to a tumor sample")
    final_normal_sample_to_tumor_sample[key] = tumor_bam_to_sample[tumor_bam]
    

sample_to_position = {}
sample_to_output = {}

header_loaded = False
line_count = 0

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
        
        for key, value in final_normal_sample_to_tumor_sample.iteritems():
            flname = open(os.path.join(args.outputdir,value +".bed"),"w")
            sample_to_output[key] = flname


    # process the rest of the lines in the file        
    else:
        line_count+= 1
        if line_count % 100000 == 0:
            print "line count == " + str(line_count)

        if header_loaded != True:
            raise NameError("Header not seen in the file")
        # check that the site is a PASS
        sp = line.strip().split("\t")
        if sp[6] != "PASS":
            continue

        # get the position
        chrom = sp[0]
        pos = sp[1]

        # check that the call is a single base (SNV) and not an indel
        for x in sp[4].split(","):
            if not x.upper() in ["A","C","G","T"]:
                continue
        
        # now check each sample; if they're a het print it to the respective file
        for key, val in sample_to_position.iteritems():
            geno_full = sp[sample_to_position[key] + 9]
            geno_each = geno_full.split(":")[0].split("/")
            if len(geno_each) != 2 or geno_each[0] == geno_each[1]:
                continue

            # ok, we have a het to output, dump it to the target file
            output_file = sample_to_output[key]
            output_file.write(chrom + "\t" + str(int(pos)-1) + "\t" + pos + "\t" + (geno_full) + "\n")


# write the output marker, pretty much a unix touch on the file
marker = open(os.path.join(args.outputdir,"complete.marker"),"w")
marker.close()
        

