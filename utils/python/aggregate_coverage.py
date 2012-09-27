# system level imports
import sys
import os
import argparse
import subprocess
import shutil
import numpy
from scipy import stats
import warnings

numpy.seterr(divide='ignore', invalid='ignore')


# create the command line arguments
parser = argparse.ArgumentParser()

# the analysis set name
parser.add_argument('-a', '--analysis_name',help="the analysis set name",required=True)

# coverage files
parser.add_argument('-t', '--tumor',help="the tumor files",required=True)
parser.add_argument('-n', '--normal',help="the normal files",required=True)

# track information
parser.add_argument('-L', '--intervals',help="the intervals bed files",required=True)
parser.add_argument('-on', '--output_normals',help="the normal output file",required=True)
parser.add_argument('-ot', '--output_tumors',help="the tumor output file",required=True)
parser.add_argument('-bf', '--output_bf',help="the bait factor",required=True)

# where to find the additional pyhton code ; in Firehose this is the libdir
parser.add_argument('-ld', '--libdir',help="the library direcrory; where to find our associated python files",required=True)
args = parser.parse_args()

# import local includes
import sys
sys.path.append(args.libdir)

from sample_handler import *
from capseg_utils import *

# ---------------------------- start processing ---------------------------
# load up the files
print "loading the normals..."
normals = load_up_tn_entries(args.normal)
print "loading the tumors..."
tumors = load_up_tn_entries(args.tumor)

# load up all the copy ratio statistics from each sample
cr_stats_normals = []
cr_stats_tumors = []

print "Getting the CR stats for each sample..."
for vals in normals.values():
    cr_stats_normals.extend(vals.cr_stats)
for vals in tumors.values():
    cr_stats_tumors.extend(vals.cr_stats)

# find the cr stat cutoff
normal_cutoff = stats.scoreatpercentile(cr_stats_normals,80)

# now tell each lane we've collected the CR stat cutoff to use
[cov.set_lanes_to_use(normal_cutoff,debugging=True) for cov in normals.values()]
[cov.set_lanes_to_use(normal_cutoff,debugging=True) for cov in tumors.values()]

print "getting the target names..."
target_file = open(args.intervals)
target_pos = []
targets = []

for line in target_file:
    sp = line.strip().split()
    targets.append(sp[3])
    target_pos.append(str(sp[0] + "\t" + sp[1] + "\t" + sp[2]))
bait_factors = []
bait_factor_names = []

# preprocess the normals
bait_factors = process_bait_factors(targets,normals,"stats_file.txt",target_pos)

# now figure out what baits to filter out
bf_cutoff_low = stats.scoreatpercentile(filter(lambda x: x > 0, bait_factors),25)

# make a list of baits to keep, and write out the bait factors
baits_to_keep = {}
bait_factors_output = open(args.output_bf,"w")
bait_factors_output_debug = open(args.output_bf + ".debug","w")
bait_factors_output.write("bait_factor\n")
bait_factors_output_debug.write("bait_factor\tkept\n")

for i in range(0,len(bait_factors)):
    if bait_factors[i] >= bf_cutoff_low: #  and bait_factors[i] <= bf_cutoff_high:
        bait_factors_output.write(targets[i] + "\t" + str(bait_factors[i]) + "\n")
        baits_to_keep[targets[i]] = True
    bait_factors_output_debug.write(targets[i] + "\t" + str(bait_factors[i]) + "\t" + str(bait_factors[i] >= bf_cutoff_low) + "\n" )

normals = load_up_tn_entries(args.normal)
# now tell each lane we've collected
[cov.set_lanes_to_use(normal_cutoff) for cov in normals.values()]

print "Processing the normal data..."
process_output(args.output_normals,bait_factors,baits_to_keep,normals,targets)
print "Processing the tumor data..."
process_output(args.output_tumors,bait_factors,baits_to_keep,tumors,targets)



