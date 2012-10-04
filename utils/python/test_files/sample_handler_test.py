import sys
sys.path.append("../")

from sample_handler import *
import unittest

class TestSampleHandler(unittest.TestCase):
    def setUp(self):
        self.coverage_file = "data/AD2031-T3.tumor.tmp.coverage.post"
        self.cr_stat_file = "data/AD2031-T3.tumor.tmp.coverage.cr.stat"
        self.column_sums = "data/AD2031-T3.tumor.tmp.coverage.col.sums"
        self.truth_data = "data/AD2031-T3.median.truth"
        self.truth_data_cut = "data/AD2031-T3.median.truth.cut"
        self.tumor_name = "AD2031-T3"

    def test_basic_reader(self):
        # create our sample handler
        sh = CoverageManager(self.coverage_file,
                             self.cr_stat_file,
                             self.tumor_name)

        # artificially set the lanes to use from the coverage manager
        sh.use_lane = [True]*len(sh.cr_stats)

        # load up the truth data
        truth = {}
        tag_order = []
        truth_reader = open(self.truth_data,"r")
        truth_reader.readline()
        for line in truth_reader:
            sp = line.strip().split("\t")
            truth[sp[0]] = float(sp[1])
            tag_order.append(sp[0])

        # now check each of the coverage calculations
        count = 0

        while True:
            # does the sample reader have data?             
            if not sh.has_data():
                break

            # check that we have the token available -- if not this indicates a problem
            self.assertIn(sh.tag,truth.keys(),msg="Unable to find token "  + sh.tag + " in the truth data")

            # now check that the values are the same
            gold_val = truth[sh.tag]
            calced = sh.get_output_value(1.0)
            
            # make sure we have the right tag
            self.assertEqual(sh.tag,tag_order[count])
            
            # make sure the median agree with the truth data
            if calced != gold_val:
                print ",".join([str(t) for t in sh.coverage])
                print sh.line
            self.assertAlmostEqual(calced, gold_val, places=7, 
                                   msg="The value [truth] " + str(gold_val) + " doesn't equal the calc'ed value " + str(calced) + " for target " + sh.tag)

            count += 1
            if count > 10000:
                break

            sh.next()

        # let the user know we looked at X number of sites
        print "Checked " + str(count) + " sites in the " + self.tumor_name

    def test_lane_elimination(self):
        # create our sample handler
        sh = CoverageManager(self.coverage_file,
                             self.cr_stat_file,
                             self.tumor_name)

        # set the lanes to use from the coverage manager
        sh.set_lanes_to_use(2.0e-6)

        # load up the truth data
        truth = {}
        tag_order = []
        truth_reader = open(self.truth_data_cut,"r")
        truth_reader.readline()
        for line in truth_reader:
            sp = line.strip().split("\t")
            truth[sp[0]] = float(sp[1])
            tag_order.append(sp[0])

        # now check each of the coverage calculations
        count = 0

        while True:
            # does the sample reader have data?             
            if not sh.has_data():
                break

            # check that we have the token available -- if not this indicates a problem
            self.assertIn(sh.tag,truth.keys(),msg="Unable to find token "  + sh.tag + " in the truth data")

            # now check that the values are the same
            gold_val = truth[sh.tag]
            calced = sh.get_output_value(1.0)
            
            # make sure we have the right tag
            self.assertEqual(sh.tag,tag_order[count])
            
            # make sure the median agree with the truth data
            if calced != gold_val:
                print ",".join([str(t) for t in sh.coverage])
                print sh.line
            self.assertAlmostEqual(calced, gold_val, places=7, 
                                   msg="The value [truth] " + str(gold_val) + " doesn't equal the calc'ed value " + str(calced) + " for target " + sh.tag)

            count += 1
            if count > 10000:
                break

            sh.next()

        # let the user know we looked at X number of sites
        print "Checked " + str(count) + " sites in the " + self.tumor_name

if __name__ == '__main__':
    unittest.main()
