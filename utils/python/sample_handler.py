import numpy
import mmap
import contextlib


'''
The CoverageManager class

This class manages coverage tracks from samples.  Each samples coverage file is loaded,
along with the CR stat.

'''
class CoverageManager:
    def __init__(self,coverage_file,cr_stat_file,tumor):
        self.cov_file = coverage_file
        self.cr_file = cr_stat_file
        self.tumor = tumor

        # load up the cr.stats and column sums for the sample
        self.cr_stats = []

        stf = open(self.cr_file)
        header = stf.readline()

        for line in stf:
            sp = line.strip().split()
            self.cr_stats.append(float(sp[1]))

        # setup the buffer of lines as a memory mapped file for speed (I hope)
        self.cov_buffer = open(self.cov_file,"r")

        self.header = self.cov_buffer.readline().strip().split("\t")

        # load the initial coverage from the first entry in the coverage file
        self.tag = "UNKNOWN"
        self.line_number = 0
        self.next()

    def get_current_tag(self):
        ''' get the current tag'''
        return self.tag

    def set_lanes_to_use(self,cr_cutoff,keep_min_percent=50,debugging=False):
        '''
        given the CR stat cutoff, figure out which lanes to drop -- the logic here is that we'd like to throw away every lane
        that has a CR stat below the threshold UNLESS we don't have enough lanes to make up the sample.  Then we're forced to
        take bad lanes, so we take enough from the best rejects to get the minimum number.

        the prereq. is that the cr_stats object has been filled in
        '''
        self.use_lane = [False]*len(self.cr_stats)
        val = []
        for i in range(0,len(self.cr_stats)):
            val.append((i,self.cr_stats[i]))
        val = sorted(val, key=lambda value: value[1])

        trues = 0
        if debugging: print "cutoff " + str(cr_cutoff)
        for i in range(0,len(self.cr_stats)):
            if self.cr_stats[i] <= cr_cutoff:
                if debugging: print "using lane " + self.header[i] + " cr stats " + str(self.cr_stats[i])
                self.use_lane[i] = True
                trues += 1
            else:
                if debugging: print "not using lane " + self.header[i]+ " cr stats " + str(self.cr_stats[i])

        ind = 0
        while trues < keep_min_percent/100.0 * len(self.use_lane):
            if debugging: print "reusing lane " + self.header[val[ind][0]]+ " cr stats " + str(self.cr_stats[val[ind][0]])
            self.use_lane[val[ind][0]] = True
            trues += 1
            ind += 1
        self.use_lane = [True]*len(self.cr_stats) # remove after testing
        print "using " + str(trues) + " of " + str(len(self.header)) + " lanes for sample " + self.cov_file

    def good_lane_count(self):
        '''return the count of the number of lanes we've kept'''
        cnt = 0
        for val in self.use_lane:
            if val:
                cnt += 1
        return cnt

    def get_coverage(self,target):
        if target != self.tag:
            raise NameError("Tags are out of sync; we expected " + self.tag + " but we were requested coverage for tag " + target)

        ''' get the raw coverage for this lane (for kept lanes)'''
        ret = []
        for i in range(0,len(self.coverage)):
            if self.use_lane[i]:
                ret.append(self.coverage[i])
        return ret

    def get_output_value(self,bait_factor,target):
        '''get the output values (the values, divided by the column sum, calibrated (divided) by the bait factor'''
        if target != self.tag:
            raise NameError("Tags are out of sync; we expected " + self.tag + " but we were requested coverage for tag " + target)

        ret = []
        for i in range(0,len(self.coverage)):
            if self.use_lane[i]:
                ret.append(self.coverage[i]/bait_factor)

        # check that we have at least one value
        if len(ret) <= 0:
            return None
        elif len(ret) == 1:
            return ret[0]
        else:
            return numpy.median(ret)

    def next(self):
        self.line = self.cov_buffer.readline()
        if not self.line:
            self.tag = None
            self.coverage = []
        sp = self.line.strip("\n").split("\t")
        self.tag = sp[0]

        self.coverage = [float(s) for s in sp[1:(len(sp))]]
        self.line_number += 1
        if len(self.coverage) != len(self.header):
            self.tag = None
            self.coverage = []
            #raise NameError("Sizes don't match; line \"" + self.line + "\" with size " + str(len(self.line)) + " (line " + str(self.line_number) + " ) doesn't have " + str(len(self.header)) + " tokens  for file " + self.cov_file)

    def has_data(self):
        if self.tag != None:
            return True
        return False
