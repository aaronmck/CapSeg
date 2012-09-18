'''
this scripts takes a firehose workspace and individual set, and creates a CapSeg run setup for the
samples
'''
import os
import argparse
import sys
import shutil
import getpass

from FirehoseConnection import *

def createFirehoseRunScript(base_dir,
                            capseg_base_dir,
                            results_directory,
                            temperary_directory,
                            sample_file,
                            run_script,
                            reference_pack,
                            analysis_set,
                            short_queue="hour",
                            long_queue="week",
                            interval_per_split=400):
    ''' create a firehose run file, from a variety of settings we have, and dump it to the run_script location'''
    name = analysis_set[0]
    if len(analysis_set) > 1:
        name += "_plus_" + str(len(analysis_set)) + "_others"
    rs = open(run_script,"w")
    line_ending = " \\" + "\n"
    rs.write("java -Djava.io.tmpdir=" + os.path.abspath(temperary_directory) + line_ending)
    rs.write("-Xmx8g \\" + "\n")
    rs.write("-jar " + os.path.join(capseg_base_dir,"utils/Queue.jar") + line_ending)
    rs.write("-S " + os.path.join(capseg_base_dir,"utils/CapSeg.scala") + line_ending)
    rs.write("-I " + sample_file + line_ending)
    rs.write("-L " + reference_pack.interval_list + line_ending)
    rs.write("-T " + reference_pack.bed_file + line_ending)
    rs.write("--tmpdir " + temperary_directory + line_ending)
    rs.write("-R " + reference_pack.reference + line_ending)
    rs.write("-bsub" + line_ending)
    rs.write("-run" + line_ending)
    rs.write("-pl" + line_ending)
    rs.write("-build " + reference_pack.build + line_ending)
    rs.write("-o " + results_directory + line_ending)
    rs.write("-jobQueue " + short_queue + line_ending)
    rs.write("-jobProject CAPSEG_" + name + line_ending)
    rs.write("-memLimit 2" + line_ending)
    rs.write("--log_to_file " + os.path.join(results_directory,"logFile.txt") + line_ending)
    rs.write("-iPerSplit " + str(interval_per_split) + line_ending)
    rs.write("--longqueue " + long_queue + line_ending)
    rs.write("--baitcsv " + reference_pack.csv_file + line_ending)
    rs.write("-retry 2" + line_ending)
    rs.write("-asn " + name + line_ending)
    rs.write("-lib " + capseg_base_dir + line_ending)
    rs.close()

class ReferencePack:
    ''' a class that stores all of the reference specific information'''
    def __init__(self,build="hg19"):
        self.build = build
        if build == "hg19":
            self.reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
            self.csv_file = "/xchip/cga2/aaron/copy_number/data/targets/hg19/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.csv"
            self.bed_file = "/xchip/cga2/aaron/copy_number/data/targets/hg19/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.bed"
            self.interval_list = "/xchip/cga2/aaron/copy_number/data/targets/hg19/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list"
        elif build == "hg18":
            self.reference = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"
            self.csv_file = "/xchip/cga2/aaron/copy_number/data/targets/hg18/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly18.targets.csv"
            self.bed_file = "/xchip/cga2/aaron/copy_number/data/targets/hg18/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly18.targets.bed"
            self.interval_list = "/xchip/cga2/aaron/copy_number/data/targets/hg18/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly18.targets.interval_list"
        elif build == "mouse":
            self.reference = "/seq/references/Mus_musculus_assembly9/v1/Mus_musculus_assembly9.fasta"
            self.csv_file = "/xchip/cga2/aaron/copy_number/data/targets/mm9/mouse_whole_exome_agilent_v1.targets.csv"
            self.bed_file = "/xchip/cga2/aaron/copy_number/data/targets/mm9/mouse_whole_exome_agilent_v1.targets.bed"
            self.interval_list = "xchip/cga2/aaron/copy_number/data/targets/mm9/mouse_whole_exome_agilent_v1.targets.interval_list"
        else:
            raise NameError("Unknown build type: " + build + ", please specify a known build type (hg19, hg18, mouse)")

# not really required, but who knows who will source this
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Gather all of the individuals, their normal names, and their cleaned normal bam files into one list.')
    parser.add_argument('--output_location',help="the output location",required=True)
    parser.add_argument('--workspace',help="the workspace",required=True)
    parser.add_argument('--individual_set',help="the individual set; this can be specified multiple times to include multiple sets together",required=True,action='append')
    parser.add_argument('--temp_space',help="where we can put temperary files (these can be large, hptmp at the Broad is recommended)",required=True)
    parser.add_argument('--capseg_base_dir',help="the base directory for CapSeg (the directory should have an R subdirectory, a utils subdir, etc)",required=True)
    parser.add_argument('--master_table_name',help="the master table name to make in the output directory (not required, defaults to sample.table.txt)",default="sample.table.txt")
    parser.add_argument('--long_queue',help="where to send long running jobs to (default week)",default="week")
    parser.add_argument('--short_queue',help="where to send short running jobs to (default hour).  You can try hour first, but if jobs fail set this to week",default="hour")
    parser.add_argument("--build",help="the build to use.  Either hg18, hg19, or mouse currently",required=True)
    parser.add_argument('--overwrite_existing_files',help="overwrite any existing files",default=False,action="store_true")
    parser.add_argument('--verbose',help="dump more information about the sample processing done with Firehose",default=False,action="store_true")
    args = parser.parse_args()

    # ------------ setup a bunch of the files, and check some pre-reqs ------------

    # the master table we're going to write to
    master_table = os.path.join(args.output_location,args.master_table_name)

    if not os.path.exists(args.output_location):
        raise NameError("Unable to find the output path " + args.output_location + ", make sure it exists")

    if os.path.exists(master_table) and not args.overwrite_existing_files:
        raise NameError("The output master_table_name already exists, specify --overwrite_existing_files to override this")

    # create the needed directories if they don't exist, and check the tmp space exists
    tmp_dir = args.temp_space
    results_dir = os.path.join(args.output_location,"results")
    if not os.path.exists(tmp_dir):
        raise NameError("unable to access the temperary space: " + tmp_dir)
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    # ------------ now run, getting the firehose samples, and creating the output ------------

    # now get their username and password - we don't want to store this in their bash history so we have to use the getpass package
    user = getpass.getuser()
    password = getpass.getpass("Firehose password: ")

    # setup a connection
    connection = FirehoseConnection(user,password)
    bam_files = connection.mappingToOutputFile(args.workspace,args.individual_set,master_table)

    # finally, make the script to run the whole thing
    createFirehoseRunScript(args.output_location,
                            args.capseg_base_dir,
                            results_dir,
                            tmp_dir,
                            master_table,
                            "run.sh",
                            ReferencePack(args.build),
                            args.individual_set,
                            short_queue=args.short_queue,
                            long_queue=args.long_queue)
