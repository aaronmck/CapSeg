# --------------------------------------------------------------
# load up the data from the tumor or normal CSV file
def load_up_tn_entries(fl):
    tn_files = open(args.tumor)
    hdr = tn_files.readline()
    ret = dict()
    for line in tumor_files:
        sp = line.strip().split("\t")
        print sp[0] + ",",
        ret[sp[0]]=CoverageManager(sp[1],sp[2],sp[3],True)
    return ret

# --------------------------------------------------------------
# process the bait factor values
def process_bait_factors(targets,coverage_managers):
    bait_factors = []
    for target in targets:
        coverage = []
        # for the given target, ask the coverage manager for the coverage for this sample
        for sample,cov_manager in normals.iteritems():
            if cov_manager.get_current_tag() == target:
                coverage.extend(cov_manager.get_coverage())
                cov_manager.next()
            else:
                coverage.extend([0]*cov_manager.good_lane_count())

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
        if processed % 5000 == 0:
            print "Preprocessed " + str(processed) + " sites"
    stats_file.close()
    return bait_factors

# --------------------------------------------------------------
# process a file, given the bait factor
def process_output(output_file,bait_factors,coverage_managers):
    output = open(output_file,"w")
    output.write("\t".join(tumors.keys()) + "\n")
    processed = 0
    for i in range(0,len(targets)):
        output_values = []
        bf = bait_factors[i]
        for sample,cov_manager in coverage_managers.iteritems():
            if cov_manager.get_current_tag() == target:
                output_values.append(cov_manager.get_output_value(bf))
                cov_manager.next()
            else:
                output_values.append(0)

        target = targets[i]
        if not baits_to_keep.has_key(target):
            continue

        processed += 1
        if processed % 10000 == 0:
            print "Processed " + str(processed) + " tumor sites"
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
