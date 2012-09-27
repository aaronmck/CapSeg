from sample_handler import *
# system level imports
import sys
import os
import argparse
import subprocess
import shutil
import numpy
from scipy import stats
import warnings

# --------------------------------------------------------------
# load up the data from the tumor or normal CSV file
def load_up_tn_entries(fl):
    tn_files = open(fl)
    hdr = tn_files.readline()
    ret = dict()
    for line in tn_files:
        sp = line.strip().split("\t")
        print sp[0] + ",",
        ret[sp[0]]=CoverageManager(sp[1],sp[2],sp[3],True)
    print "done loading file: " + fl
    return ret

# --------------------------------------------------------------
# process the bait factor values
def process_bait_factors(targets,coverage_managers,stats_filename,target_pos):
    stats_file = open(stats_filename,"w")
    stats_file.write("target\tbf\tentropy\tvar\tcontig\tstart\tstop\n")

    bait_factors = []
    index = 0
    processed = 0
    for target in targets:
        coverage = []
        # for the given target, ask the coverage manager for the coverage for this sample
        for sample,cov_manager in coverage_managers.iteritems():
            if cov_manager.get_current_tag() == target:
                coverage.extend(cov_manager.get_coverage())
                cov_manager.next()
            else:
                coverage.extend([0]*cov_manager.good_lane_count())
                print "missing coverage for target " + target + " in sample " + sample + " at index " + str(index)

        # calculate a couple of statistics on each bait
        bf = numpy.median(coverage)
        entropy = stats.distributions.entropy(coverage)
        sm = numpy.sum(coverage)
        variance = numpy.var([s/sm for s in coverage])

        # write to the stat file
        stats_file.write(target + "\t" + str(bf) + "\t" + str(entropy) + "\t" + str(variance) + "\t" + target_pos[index] + "\n")
        index += 1
        bait_factors.append(bf)
        processed += 1
        if processed % 10000 == 0:
            print "Preprocessed " + str(processed) + " sites"
    stats_file.close()
    return bait_factors

# --------------------------------------------------------------
# process a file, given the bait factor
def process_output(output_file,bait_factors,baits_to_keep,coverage_managers,targets):
    output = open(output_file,"w")
    output.write("\t".join(coverage_managers.keys()) + "\n")
    processed = 0
    for i in range(0,len(targets)):
        target = targets[i]
        output_values = []
        bf = bait_factors[i]

        for sample,cov_manager in coverage_managers.iteritems():
            if cov_manager.get_current_tag() == target:
                output_values.append(cov_manager.get_output_value(bf,target))
                cov_manager.next()
            else:
                output_values.append(0)

        if not baits_to_keep.has_key(target):
            continue

        processed += 1
        if processed % 10000 == 0:
            print "Processed " + str(processed) + " sites"
        output.write(target + "\t")
        output.write("\t".join([str(s) for s in output_values]) + "\n")

# --------------------------------------------------------------
# write out the bait factors
def write_out_bait_factors(bf_out_file,baits_to_keep,targets):
    bf_out = open(bf_out_file,"w")
    bf_out.write("bait\tbait.factor\n")
    for i in range(0,len(targets)):
        target = targets[i]
        if not baits_to_keep.has_key(target):
            continue
        bf_out.write(target + "\t" + str(bait_factors[i]) + "\n")
    bf_out.close()
