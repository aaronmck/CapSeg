import numpy
'''
The CoverageManager class

This class manages coverage tracks from samples.  Each samples coverage file is loaded,
along with the CR stat.

'''
class CoverageManager:
    def __init__(self,coverage_file,cr_stat_file,column_sums,tumor):
        self.cov_file = coverage_file
        self.cr_file = cr_stat_file
        self.tumor = tumor

        # load up the cr.stats and column sums for the sample
        self.cr_stats = []

        stf = open(self.cr_file)
        header = stf.readline()
        print self.cr_file
        for line in stf:
            sp = line.strip().split()
            self.cr_stats.append(float(sp[1]))

        # setup the buffer of lines
        self.cov_buffer = open(self.cov_file)
        self.header = self.cov_buffer.readline().strip().split("\t")

        # load the initial coverage from the first entry in the coverage file
        self.tag = "UNKNOWN"
        self.next()

    def get_current_tag(self):
        ''' get the current tag'''
        return self.tag

    def set_lanes_to_use(self,cr_cutoff,keep_min_percent=50,debugging=False):
        '''given the CR stat cutoff, figure out wgat lanes to drop'''
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
        print "using " + str(trues) + " of " + str(len(self.header)) + " lanes for sample " + self.cov_file
        self.use_lane = [True]*len(self.cr_stats)

    def good_lane_count(self):
        '''return the count of the number of lanes we've kept'''
        cnt = 0
        for val in self.use_lane:
            if val:
                cnt += 1
        return cnt

    def get_coverage(self):
        ''' get the raw coverage for this lane (for kept lanes)'''
        ret = []
        for i in range(0,len(self.coverage)):
            if self.use_lane[i]:
                ret.append(self.coverage[i])
        return ret

    def get_output_value(self,bait_factor,target):
        '''get the output values (the values, divided by the column sum, calibrated (divided) by the bait factor'''
        ret = []
        # print "bait factor " + str(bait_factor)
        for i in range(0,len(self.coverage)):
            if self.use_lane[i]:
                ret.append(self.coverage[i]/bait_factor)
        med = ret[0]
        # the median code doesn't like lengths of [0-1) so we might as well save time and check for [0-1]
        if len(ret) > 1:
            med = numpy.median(ret)
            #print target
            #print self.tag
            #print "tot = " + "\t".join([str(r) for r in self.coverage])
            #print "sample " + self.cov_file + "\t" + "\t".join([str(r) for r in ret]) + " med " + str(med)
        return med

    def next(self):
        self.line = self.cov_buffer.readline()
        sp = self.line.strip().split("\t")
        self.tag = sp[0]
        self.coverage = [float(s) for s in sp[1:(len(sp))]]

