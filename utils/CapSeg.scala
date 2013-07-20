import scala.util._
import java.lang.Integer
import org.broadinstitute.sting.commandline.ArgumentSource
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.function.scattergather.{GatherFunction, CloneFunction, ScatterFunction}
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import collection.JavaConversions._
import org.broadinstitute.sting.utils.interval.IntervalUtils
import org.broadinstitute.sting.commandline.{Argument, Output, Input}
import org.broadinstitute.sting.queue.function.CommandLineFunction
import java.io.File
import scala.io.Source
import scala.collection.JavaConversions._
import java.util.HashMap
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.utils.interval._
import java.io._

// --------------------------------------------------------------------------------------------------------------------------------
// CapSeg: the capture based segmentation
//
// --------------------------------------------------------------------------------------------------------------------------------
class CapSeg extends QScript {
  qscript =>
  // --------------------------------------------------------
  // required arguments
  // --------------------------------------------------------
  @Input(doc = "a file containing a sample column, a tumor column, and a normal column", shortName = "I")
  var samples: File = _

  @Argument(doc = "path to the reference", shortName = "R")
  var reference: File = _

  @Input(doc = "The path to the csv list of the targets", shortName = "bed_csv")
  var baitCSV: File = _

  @Argument(doc = "The BED file containing the intervals for each of the targets captured in WES", shortName = "T")
  var targets: File = _

  @Argument(doc = "the genome build: either hg18, hg19, or b36, or mm9", shortName = "build")
  var build: String = _

  @Input(doc = "The interval file containing the intervals for each of the targets captured in WES", shortName = "L")
  var intervals: File = _

  @Input(doc = "Where to store the output of the run", shortName = "o")
  var outputDir: File = _

  @Argument(doc = "our temporary space; where to store temporary files", required = true)
  var tmpDir = new File("/broad/shptmp/aaron/tmpbed/")

  @Argument(doc = "the library directory", shortName="lib", required = true)
  var libraryDir: File = _

  @Argument(doc = "what to call this analysis set", shortName = "asn", required = true)
  var analysisSet: String = _

  @Input(doc = "dbSNP file; used to order the SNPs in the allele coverage pulldown list ", shortName = "dbsnp", required = true)
  var dbsnp: File = _

  // ---------------------------------- optional arguments ----------------------------------
  @Argument(doc = "what long queue to submit jobs to", required = false)
  var longQueue = "week"

  @Argument(doc = "what short queue to submit jobs to", required = false)
  var shortQueue = "hour"

  @Argument(doc = "should we use historical planes for data?", shortName = "uhd", required = false)
  var useHistoricalData: String = _

  @Input(doc = "If we need to load the prototypical coverage for normals (we're running without normals), what file should we use?", shortName = "pn", required=false)
  var protoNormalCoverage: File = new File("/xchip/cga2/aaron/static/prototypical_normal/gbm_prototype.tmp.coverage")

  @Argument(doc = "do we want to collect coverage at the lane level (true) or at the sample level (default)", shortName = "pl", required = false)
  var perLane = false // the logic has to be reversed a little to work with the arg system

  @Input(doc = "tangent normalization input database location", shortName = "tnd", required = false)
  var tangentLocation = new File("/xchip/cga2/aaron/static/normal_subspaces/")

  @Input(doc = "tangent normalization output database location", shortName = "tndo", required = false)
  var tangentOutputLocation = new File("/xchip/cga2/aaron/static/normal_subspaces/")

  // ------------- segmentation parameters; used to when segmenting the signal for each sample ---------------
  @Argument(doc = "the alpha value to use for the cutoff in segmentation", shortName = "sav", required = false)
  var segmentAlpha = "0.001"

  @Argument(doc = "the undo splits method", shortName = "ssm", required = false)
  var segmentSplits = "sdundo"

  @Argument(doc = "the segment SD setting", shortName = "ssd", required = false)
  var segmentSD = "1.5"

  @Argument(doc = "the threshold to merge segments: segments smaller than this size with be merged with the appropriate neighbor", shortName = "", required = false)
  var segMergeThresh = 1.5

  @Argument(doc = "the probability threshold for Bayesian segment merging", shortName = "msc", required = false)
  var minSegCount = 10

  // the command line traits -- defaults for our pipeline
  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = reference
    this.memoryLimit = 6
  }

  // where to put the queue log directory
  val queueLogDir: String = ".qlog/"  // Gracefully hide Queue's output

  // ----------------------------------------------------------------------------------------------------------------
  // The main script
  // ----------------------------------------------------------------------------------------------------------------
  def script = {
    // Due to Firehose parameters we have to take in string represented boolean parameters and convert to the Boolean type inside of this script
    var useHistData = false
    if (useHistoricalData != null && useHistoricalData.toUpperCase().equals("TRUE"))
      useHistData = true

    // load up the samples
    val sampleObj = new TumorNormalFile(samples)

    // output the tumor and normal samples -- used by the sample to lane walker
    val tumorBamFileList = new File(outputDir.getAbsolutePath() + "/" + analysisSet + ".tumors.list")
    val normalBamFileList = new File(outputDir.getAbsolutePath() + "/" + analysisSet + ".normals.list")
    val tumorBams = sampleObj.writeTumorsFile(tumorBamFileList)
    val normalBams = sampleObj.writeNormalsFile(normalBamFileList)

    // setup our coverage files; these contain a line for each sample, with its coverage, cr.stat, and column sum files
    val normalCoverageFiles = new File(outputDir.getAbsolutePath() + "/normalCoverageFiles.txt")
    val tumorCoverageFiles = new File(outputDir.getAbsolutePath() + "/tumorCoverageFiles.txt")

    // pull down the coverage at all of the targets
    var normals = List.empty[String]
    if (sampleObj.getNormalMap().size == 0 && protoNormalCoverage != null)
      normals = prototypicalNormalPulldownCoverage(protoNormalCoverage,normalCoverageFiles)
    else if (sampleObj.getNormalMap().size == 0)
      throw new IllegalArgumentException("When there are no normals, make sure to set the prototypical normal parameter (see the help with -h)")
    else
      normals = pulldownCoverage(sampleObj.getNormalMap(),normalCoverageFiles,"normal")
    val tumors = pulldownCoverage(sampleObj.getTumorMap(),tumorCoverageFiles,"tumor")

    // make the build specific interval file
    var sampleIntervals = new File(tmpDir.getAbsolutePath() + "/sample.interval_list")

    var out = new java.io.FileWriter(sampleIntervals)
    if (build.equals("hg18") || build.equals("mm9")) { out.write("chr1:1\n")} else { out.write("1:1\n") }
    out.close()

    //  we have four files to create -- sample to lane, and bam to sample files for each tumor and normal
    val tumorSampleToPGFile = new File(outputDir.getAbsolutePath() + "/tumorSampleInformation.txt")
    val tumorBamToSampleFile = new File(outputDir.getAbsolutePath() + "/tumorBamToSample.txt")
    val normalSampleToPGFile = new File(outputDir.getAbsolutePath() + "/normalSampleInformation.txt")
    val normalBamToSampleFile = new File(outputDir.getAbsolutePath() + "/normalBamToSample.txt")

    // add the sample to lane walker - get the sample to file and bam to sample tables for both tumor and the normal
    sampleToLane(tumorBamFileList, tumorSampleToPGFile, tumorBamToSampleFile, sampleIntervals)
    if (sampleObj.getNormalMap().size > 0)
      sampleToLane(normalBamFileList, normalSampleToPGFile, normalBamToSampleFile, sampleIntervals)

    // setup our coverage files; these contain a line for each sample, with its coverage, cr.stat, and column sum files
    val normalSampleFile = new File(outputDir.getAbsolutePath() + "/normalSampleCoverage.txt")
    val tumorSampleFile = new File(outputDir.getAbsolutePath() + "/tumorSampleCoverage.txt")
    val finalNormalMatrixBF = new File(outputDir.getAbsolutePath() + "/baitFactors.txt")
    val signalFile = new File(outputDir.getAbsolutePath() + "/signalFiles.txt")

    // merge all the output into a single file
    add(new MergeCoverage(analysisSet,
                          tumorCoverageFiles,
                          normalCoverageFiles,
                          targets,
                          tumorSampleFile,
                          normalSampleFile,
                          libraryDir,
                          finalNormalMatrixBF,
                          tumors,
                          normals,
                          longQueue))

    // the jobs are tied by inputs and outputs -- this step is actually after the post process step, but we do it here since the signal file
    // names are made sample-by-sample here
    var signal = List.empty[File]
    var sampleToSignal = new java.io.FileWriter(signalFile)
    var firehoseInport = new java.io.FileWriter(analysisSet + ".segments.tsv")
    sampleToSignal.write("sample\tsignal.file\n")
    firehoseInport.write("individual_id\tcapseg_segmentation_file\n")

    // where to put the allele balance information
    val alleleFreq = new File(outputDir + "/alleleFreq/")
    alleleFreq.mkdir();

    var srcDir = new File(libraryDir.getAbsolutePath() + "/R/allelic/")
    var outAllelicDir = new File(outputDir + "/allelicCapSeg/")
    var alleleOutput = List[File]()

    // for each of the tumor samples, create segments, and for those samples with associated normal tissue,
    // run allelic CapSeg
    var normalMap = sampleObj.getNormalMap
    var vcfMap = sampleObj.getVCFMap
    for ((sample, tumor) <- sampleObj.getTumorMap) {
      val signalFile = new File(outputDir.getAbsolutePath() + "/signal/" + sample + ".tsv")
      val segmentFile = new File(outputDir.getAbsolutePath() + "/segments/" + sample + ".seg.txt")

      val output = alleleFreq + "/" + sample + ".acov"
      alleleOutput ::= output

      add(new SegmentSample(libraryDir,
                            new File(tumor).getName,
                            tumorBamToSampleFile,
                            segmentFile,
                            signalFile,
                            segmentAlpha,
                            segmentSplits,
                            segmentSD))

      // now check to ensure that this tumor has a normal sample associated with it
      if (sampleObj.checkForTumorNormalBamsAndVCF(sample)) {
        alleleBalance(tumor, normalMap(sample), vcfMap(sample), dbsnp, output, reference)
        add(new AllelicCapSeg(signalFile, outAllelicDir, segmentFile, output, tumor, srcDir, libraryDir, tumorBamToSampleFile, segMergeThresh, minSegCount, queueLogDir))
      }

      signal ::= signalFile
      sampleToSignal.write(sample + "\t" + signalFile.getAbsolutePath + "\n")
      firehoseInport.write(sample + "\t" + segmentFile.getAbsolutePath + "\n")
    }
    sampleToSignal.close()
    firehoseInport.close()

    // run the final R script which puts it all together
    add(new PostProcessData(libraryDir,
                            normalSampleFile,
                            tumorSampleFile,
                            baitCSV,
                            true,
                            new File(outputDir.getAbsolutePath() + "/.cache"),
                            normalSampleToPGFile,
                            tumorSampleToPGFile,
                            outputDir,
                            tangentLocation,
                            tangentOutputLocation,
                            build,
                            analysisSet,
                            perLane,
                            finalNormalMatrixBF,
                            signal,
                            signalFile,
                            useHistData,
                            alleleOutput))
  }


  // ----------------------------------------------------------------------------------------------------------------
  // utility function section; also contains GATK walker invocations
  // ----------------------------------------------------------------------------------------------------------------

  // sample to lane walker
  def sampleToLane(bamsIn: File, sampleToLane: File, bamFileToSample: File, sampleInterval: File) = {
    val sampleTL = new SampleToLaneWalker
    sampleTL.input_file :+= bamsIn
    sampleTL.strf = sampleToLane
    sampleTL.bs = bamFileToSample
    sampleTL.reference_sequence = reference
    sampleTL.memoryLimit = 2
    sampleTL.intervals :+= sampleInterval
    add(sampleTL)
  }

  // get allelic information at the target het sites in a sample
  def alleleBalance(tumorBam: File, normalBam: File, vcf: File, dbSNP: File, output: File, ref: File) = {
    val aBal = new AlleleCountWalker
    // put in the tagged bam files
    var bams: List[File] = Nil
    bams :+= tumorBam
    bams :+= normalBam

    var bamNames: List[String] = Nil
    bamNames :+= "tumor"
    bamNames :+= "normal"

    aBal.input_file = bamNames.zip(bams).map { case (name, file) => TaggedFile(file,name) }
    aBal.DbSNP = dbSNP
    aBal.calls = vcf
    aBal.allelebalance = output
    aBal.intervals :+= vcf
    aBal.reference_sequence = reference
    add(aBal)
  }

  def pulldownCoverage(sampleMap: Map[String, File], coverageFile: String, tumorNormal: String): List[File] = {
    val sFile = new PrintWriter(new File(coverageFile))
    sFile.write("sample\tcoverage\tcr.stat\tcolumn.sums\n")
    var returnList = List[String]()
    // Get the bam files from the list
    for ((sample,bamfile) <- sampleMap) {
      var outputInitCov = tmpDir + "/" + sample + "." + tumorNormal + ".tmp.coverage"
      var outputPost = tmpDir + "/" + sample + "." + tumorNormal + ".tmp.coverage.post"
      var outputCRStat = tmpDir + "/" + sample + "." + tumorNormal + ".tmp.coverage.cr.stat"
      var outputColSums = tmpDir + "/" + sample + "." + tumorNormal + ".tmp.coverage.col.sums"
      returnList ::= new File(outputColSums)
      baitDepth(bamfile, targets, outputInitCov, false)

      add(new PostProcessBaitCoverage(outputInitCov, targets, outputPost, outputCRStat, outputColSums, libraryDir))
      sFile.write(sample + "\t" + outputPost + "\t" + outputCRStat + "\t" + outputColSums + "\n")
    }
    sFile.close()
    return(returnList)
  }

  // sometimes we don't have a normal, so we'd like to use the prototypical coverage
  def prototypicalNormalPulldownCoverage(protoLocCoverage: File, coverageFile: String): List[String] = {
    val sFile = new PrintWriter(new File(coverageFile))
    sFile.write("sample\tcoverage\tcr.stat\tcolumn.sums\n")
    var returnList = List[String]()

    var outputInitCov = protoLocCoverage.getAbsolutePath()
    var outputPost = outputInitCov + ".post"
    var outputCRStat = outputInitCov + ".cr.stat"
    var outputColSums = outputInitCov + ".col.sums"
    returnList ::= new File(outputColSums)
    sFile.write("prototypicalNormal\t" + outputPost + "\t" + outputCRStat + "\t" + outputColSums + "\n")
    sFile.close()
    return(returnList)
  }

  // sample to lane walker
  def baitDepth(bamsIn: File, bed: File, outputFile: File, perSample: Boolean) = {
    val bait = new BaitDepthWalker with CommandLineGATKArgs
    bait.input_file :+= bamsIn
    bait.bed = bed
    bait.perSample = perSample
    bait.out = outputFile
    bait.intervals :+= bed
    add(bait)
  }
}

// ----------------------------------------------------------------------------------------------------------------
// the classes - these tell the CapSeg script how to call external tools
// ----------------------------------------------------------------------------------------------------------------

// create firehose import file
class MergeCoverage(analysis: String, tumorFile: File, normalFile: File, intervalsBED: File, outputTumor: File, outputNormal: File, libraryDir: File, baitFactorFile: File, tumorFileList: List[File], normalFileList: List[File], longQueue: String) extends CommandLineFunction {
  @Input(doc = "the library directory") var libDir = libraryDir
  @Input(doc = "the tumor tab delimited file") var tFile = tumorFile
  @Input(doc = "the normal tab delimited file") var nFile = normalFile
  @Input(doc = "the intervals file") var intervals = intervalsBED
  @Input(doc = "the list of tumor files (to tell queue we're dependent on them to run)") var tumors = tumorFileList
  @Input(doc = "the list of normal files (to tell queue we're dependent on them to run)") var normals = normalFileList
  @Argument(doc = "the analysis set name") var analysisSet = analysis

  @Output(doc = "the output tumor tsv") var outT = outputTumor
  @Output(doc = "the output normal tsv") var outN = outputNormal
  @Output(doc = "the bait factor information") var bf = baitFactorFile

  // some job settings
  memoryLimit = Some(2)
  jobQueue = longQueue

  def commandLine = "python %s/utils/python/aggregate_coverage.py -a %s -t %s -n %s -L %s -ot %s -on %s -bf %s -ld %s".format(libDir.getAbsolutePath(),analysisSet,tFile.getAbsolutePath(),nFile.getAbsolutePath(),intervals.getAbsolutePath(),outT.getAbsolutePath(),outN.getAbsolutePath(),bf.getAbsolutePath(),libDir.getAbsolutePath())
}


// segment samples
class SegmentSample(libraryDir: File, bamName: String, bamToSample: File, outputFl: File, signalFl: File, alpha: String, splits: String, SD: String) extends CommandLineFunction {
  @Argument(doc = "the bam file name") var bamFile = bamName
  @Input(doc = "the bam to sample file") var bamToSampleFile = bamToSample
  @Output(doc = "the output file") var outputFile = outputFl
  @Input(doc = "the signal file") var signalFile = signalFl
  @Input(doc = "library directory") var libDir = libraryDir
  @Argument(doc = "alpha value for the segmentation") var mAlpha = alpha
  @Argument(doc = "split segment rejoining approach") var mSplits = splits
  @Argument(doc = "the SD cutoff to use") var mSD = SD

  memoryLimit = Some(4)

  def commandLine = "Rscript %s/R/individual_segmentation.R --bam.file.name %s --bam.to.sample.file %s --output.filename %s --signal.file %s --r.dir %s --alpha %s --undo.splits %s --undo.SD %s".format(libDir.getAbsolutePath(), bamFile, bamToSampleFile.getAbsolutePath(), outputFile.getAbsolutePath(), signalFile.getAbsolutePath(),libDir.getAbsolutePath(), mAlpha, mSplits, mSD)
}


// post process the data using the R script
class PostProcessData(libraryDir: File, normal_bait_coverage: File, tumor_bait_coverage: File, target_file: File,
                      useCachedData: Boolean, cachedLocation: File, sToRGFileNormals: File, sToRGFileTumors: File,
                      outputDir: File, tangentLocation: File, tangentOutputLocation: File, buildType: String, analysisSet: String,
                      byLane: Boolean, normalBaitFile: File, signalFLs: List[File], signalTSV: File, histoData: Boolean, acovFiles: List[File]) extends CommandLineFunction {
  @Input(doc = "the normal bait output file") var normalFile = normal_bait_coverage
  @Input(doc = "the tumor bait output file") var tumorFile = tumor_bait_coverage
  @Input(doc = "the target list") var targetFile = target_file
  @Input(doc = "where can we find the WES segmentation algorithm") var libDir = libraryDir
  @Input(doc = "the file containing the mapping of sample to read group (normals)") var normalSampleToReadGroupFile = sToRGFileNormals
  @Input(doc = "the file containing the mapping of sample to read group (tumors)") var tumorSampleToReadGroupFile = sToRGFileTumors
  @Input(doc = "the output location") var outputDirectory = outputDir
  @Input(doc = "tangent normalization database location") var tangentDatabase = tangentLocation
  @Input(doc = "tangent normalization database location") var tangentOutputDatabase = tangentOutputLocation
  @Input(doc = "signal tsv file") var signalTSVFile = signalTSV
  @Input(doc = "coverage files from the allele balance pulldown") var coverageFiles = acovFiles

  @Output(doc = "the signal file outputs") var signalFiles = signalFLs
  @Output(doc = "the cache location") var cacheLocation = cachedLocation

  @Argument(doc = "the build type, hg18 or hg19") var build = buildType
  @Argument(doc = "the analysis name") var analysisSetName = analysisSet
  @Argument(doc = "are we running by lane?") var byLaneData = byLane
  @Argument(doc = "the normal bait factor file") var nbf = normalBaitFile
  @Argument(doc = "should we use historical data") var uhd = histoData

  //@Output(doc="the sample information file") var sampleInfoFile = new File(outputDir.getAbsolutePath() + "/sampleInformation.txt")

  @Argument(doc = "should we use the cache file") var useCache = useCachedData
  memoryLimit = Some(16) // change me
  def commandLine = "Rscript %s/R/tangent_normalize.R --normal.lane.data %s --tumor.lane.data %s --target.list %s --use.cache %s --cache.location %s --script.dir %s --normal.sample.to.lanes.file %s --tumor.sample.to.lanes.file %s --output.location %s --tangent.database.location %s --output.tangent.database %s --build %s --analysis.set.name %s --bylane %s --bait.factor %s --signal.files %s --histo.data %s".format(libDir.getAbsolutePath(), normalFile.getAbsolutePath(), tumorFile.getAbsolutePath(), targetFile.getAbsolutePath(), useCache, cacheLocation.getAbsolutePath(), libDir.getAbsolutePath(), normalSampleToReadGroupFile.getAbsolutePath(), tumorSampleToReadGroupFile.getAbsolutePath(), outputDirectory.getAbsolutePath(), tangentDatabase.getAbsolutePath(), tangentOutputDatabase.getAbsolutePath(), build, analysisSetName,byLaneData, nbf.getAbsolutePath(), signalTSVFile.getAbsolutePath(), uhd)
}

// correct the output files for any dups lines and extra headers
class PostProcessBaitCoverage(inputFile: File, targets: File, outputFile: File, outputCRSTAT: File, outputColSums: File, utilityLoc: File) extends CommandLineFunction {
  @Input(doc = "csv file outputted from a previous step") var inputCSV = inputFile
  @Input(doc = "our targets, as a bed file") var targetFile = targets

  @Output(doc = "the final csv output file") var outputCSV = outputFile
  @Output(doc = "the output cr stat file") var outputCR = outputCRSTAT
  @Output(doc = "the output column sums") var outputCS = outputColSums
  @Argument(doc = "where to find the R scripts") var loc = utilityLoc
  memoryLimit = Some(4)

  def commandLine = "Rscript %s/R/proportional_coverage.R --input.file %s --intervals %s --output.file %s --output.cr.stat.file %s --output.column.sums %s".format(loc.getAbsolutePath(),inputCSV.getAbsolutePath(),targetFile.getAbsolutePath(),outputCSV.getAbsolutePath(),outputCR.getAbsolutePath(),outputCS.getAbsolutePath())
}

// get allelic capseg data
class AllelicCapSeg(probeFL: File, outputDir: File, segmentFile: File, covFile: File, bamName: File, codeDir: File, baseScript: File, bamToSmp: File, segMerge: Double, minSeg: Integer, queueLogDir: File) extends CommandLineFunction {
  @Input(doc = "the probe file") var probe = probeFL
  @Argument(doc = "code directory; where to find the rest of the R code") var code = codeDir
  @Argument(doc = "segment merge threshold") var segMergeThresh = segMerge
  @Argument(doc = "minimum segment size") var minSegCount = minSeg
  @Input(doc = "output directory") var out = outputDir
  @Input(doc = "the segment file") var seg = segmentFile
  @Input(doc = "the coverage file") var cov = covFile
  @Input(doc = "the bam name") var bam = bamName
  @Input(doc = "the bam to sample file") var bamToSample = bamToSmp
  @Argument(doc = "where to find the base directory") var loc = baseScript
  memoryLimit = Some(2)
  this.analysisName = queueLogDir + "/" + bamName.getName() + ".allelicResults"
  this.jobName = queueLogDir + "/" + bamName.getName() + ".allelicResults"
  def commandLine = "Rscript %s/R/allelic_capseg.R --output.dir %s --probe.file %s --segment.file %s --coverage.file %s --bam.name %s --source.directory %s --bam.to.sample %s --seg.merge.thresh %s --min.seg.size %s".format(loc.getAbsolutePath(),out.getAbsolutePath(),probe.getAbsolutePath(),seg.getAbsolutePath(),cov.getAbsolutePath(),bam,code.getAbsolutePath(),bamToSample.getAbsolutePath(),segMergeThresh,minSegCount);
}

// ------------------------------------------------------------------
// -------------------- utility classes -----------------------------
// ------------------------------------------------------------------
// the base class for reading in a delimited file
class ReadDelimitedFile(inputFile: String, delimiter: String, header: Boolean = true, headerCheckString: String = "") {
  var columns = List.empty[String]
  var rows = List.empty[String]
  val delim = delimiter
  val source = Source.fromFile(new File(inputFile)).getLines
  var headerLine:String = null
  if (header && source.hasNext) {
    headerLine = source.next()
    if (!(headerLine.trim == headerCheckString.trim)) {
      throw new UserException.BadInput("Your file " + inputFile + " doesn't have the expected header: " + headerCheckString.trim + ", instead we saw " + headerLine.trim);
    }
  }
  while (source.hasNext) {
    var line = source.next()
    rows = line :: rows
  }
  //print(rows.length)
}

// ------------------------------------------------------------------
// --------------------    Utilities functions   --------------------
// ------------------------------------------------------------------
// delete a file; this waits for the trigger file before running
class DeleteMeFunction extends CommandLineFunction {
  @Argument(doc = "The file to be deleted") var me: File = _
  @Input(doc = "The file which must exist before we are allowed to delete") var trigger: File = _

  def commandLine = "rm -f %s".format(me.getAbsolutePath)
}

// delete a file; this waits for the trigger file before running
class DoubleTriggerDeleteMeFunction extends CommandLineFunction {
  @Argument(doc = "The file to be deleted") var me: File = _
  @Input(doc = "The file which must exist before we are allowed to delete") var trigger: File = _
  @Input(doc = "The second file which must exist before we are allowed to delete") var trigger2: File = _

  def commandLine = "rm -f %s".format(me.getAbsolutePath)
}

// our sample files - get tumor normal pairs
class TumorNormalFile(inputFile: String) extends ReadDelimitedFile(inputFile, "\t", true, "sample\ttumor_bam\tnormal_bam\tvcf_file\n") {
  var tBamFiles = Map.empty[String, File]
  var nBamFiles = Map.empty[String, File]
  var vcfFiles = Map.empty[String, File]
  for (line <- rows) {
    var sp = line.split(delim)
    //println(sp.length)
    if (sp.length >= 4 && sp(0) != null && sp(1) != null && sp(2) != null) {
      if (!sp(1).equals("NA") && !sp(1).equals("")) {
        tBamFiles += sp(0) -> new File(sp(1))
      }
      if (!sp(2).equals("NA") && !sp(2).equals("")) {
        nBamFiles += sp(0) -> new File(sp(2))
      }
      if (!sp(3).equals("NA") && !sp(3).equals("")) {
        vcfFiles += sp(0) -> new File(sp(3))
      }

    } else {
      throw new IllegalArgumentException("Unable to properly split line " + line)
    }
    // println("done")
  }

  // get the tumors
  def getTumors(): List[File] = {
    return (tBamFiles.valuesIterator.toList)
  }

  // get the tumors as a map
  def getTumorMap(): Map[String, File] = {
    return (tBamFiles)
  }

  // get the normal list
  def getNormals(): List[File] = {
    return (nBamFiles.valuesIterator.toList)
  }

  // get the normal samples as a map to name from file
  def getNormalMap(): Map[String, File] = {
    return (nBamFiles)
  }

  // get the vcfs as a map to name from file
  def getVCFMap(): Map[String, File] = {
    return (vcfFiles)
  }

  def printToFile(f: java.io.File)(op: java.io.PrintWriter => Unit) {
    val p = new java.io.PrintWriter(f)
    try { op(p) } finally { p.close() }
  }
  def writeTumorsFile(filename: String) {
    printToFile(new File(filename))(p => {
      getTumors().foreach(p.println)
    })
  }

  def writeNormalsFile(filename: String) {
    printToFile(new File(filename))(p => {
      getNormals().foreach(p.println)
    })
  }

  def checkForTumorNormalBamsAndVCF(tumorSample: String): Boolean = {
    // check that we have the sample.  If we do, check that we have
    // a tumor and normal bam file and a VCF file. If so return true, otherwise
    // false
    if ((getTumorMap() contains tumorSample) && (getNormalMap() contains tumorSample) && (vcfFiles contains tumorSample)) {
      return true
    }
    return false
  }
}

