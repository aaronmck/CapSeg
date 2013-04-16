import scala.util._
import org.broadinstitute.sting.commandline.ArgumentSource
import org.broadinstitute.sting.gatk.DownsampleType
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

  @Input(doc = "the VCF file with calls for the normal tissue sample in our set", shortName = "vcf", required = true)
  var vcfFile: File = _

  @Input(doc = "dbSNP file; used to order the SNPs in the allele coverage pulldown list ", shortName = "dbsnp", required = true)
  var dbsnp: File = _

  // ---------------------------------- optional arguments ----------------------------------
  @Argument(doc = "intervals per split", shortName = "iPerSplit", required = false)
  var splitSize = 150

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

  // the command line traits -- defaults for our pipeline
  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = reference
    this.intervals = List(qscript.intervals)
    this.memoryLimit = 4
  }

  // ----------------------------------------------------------------------------------------------------------------
  // The main script
  // ----------------------------------------------------------------------------------------------------------------
  def script = {
    // Firehose - we have to take in string boolean parameters and convert to the Boolean type inside of this script
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

    // create a list of the output files, one bed for each output
    val depFile = new File(alleleFreq + "/complete.marker")
    val sexFile = new File(alleleFreq + "/sex.calls")
    add(new VCFToBED(vcfFile, alleleFreq, normalBamToSampleFile, tumorBamToSampleFile, samples, libraryDir, depFile, sexFile))

    // for each tumor convert the associated normal sample BED file from the above step with the bam into an allele balance file
    var acovFiles = addAlleleBalance(sampleObj.getTumorMap(),dbsnp, alleleFreq, depFile)

    for ((sample, tumor) <- sampleObj.getTumorMap) {
      val signalFile = new File(outputDir.getAbsolutePath() + "/signal/" + sample + ".tsv")
      val segmentFile = new File(outputDir.getAbsolutePath() + "/segments/" + sample + ".seg.txt")
      add(new SegmentSample(libraryDir,
                            new File(tumor).getName,
                            tumorBamToSampleFile,
                            segmentFile,
                            signalFile))
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
                            acovFiles,
                            sexFile))
  }

  // sample to lane walker
  def sampleToLane(bamsIn: File, sampleToLane: File, bamFileToSample: File, sampleInterval: File) = {
    val sampleTL = new SampleToLane
    sampleTL.input_file :+= bamsIn
    sampleTL.strf = sampleToLane
    sampleTL.bs = bamFileToSample
    sampleTL.reference_sequence = reference
    sampleTL.memoryLimit = 2
    sampleTL.intervals :+= sampleInterval
    add(sampleTL)
  }

  // get allelic information at the target het sites in a sample
  def alleleBalance(bam: File, bed: File, dbSNP: File, output: File, dep: File) = {
    val aBal = new AlleleCount
    aBal.input_file :+= bam
    aBal.DbSNP = dbSNP
    aBal.calls = bed
    aBal.output = output
    aBal.dependency = dep
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

  // for each tumor sample, pull down the appropriate allele balance
  def addAlleleBalance(tumorMap: Map[String, File], dbsnp: File, alleleFreq: File, depFile: File): List[File] = {
    var returnList = List[File]()
    for ((sample,bamfile) <- tumorMap) {
      println("sample = " + sample + " bam = " + bamfile)
      val output = alleleFreq + "/" + sample + ".acov"
      // now figure out what the normal BED file is called
      val bedfile = alleleFreq + "/" + sample + ".bed"
      alleleBalance(bamfile, bedfile, dbsnp, output, depFile)
      returnList ::= output
    }
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
    val bait = new BaitDepth with CommandLineGATKArgs
    bait.input_file :+= bamsIn
    bait.bed = bed
    bait.perSample = perSample
    bait.out = outputFile
    bait.isr = IntervalSetRule.INTERSECTION
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
class SegmentSample(libraryDir: File, bamName: String, bamToSample: File, outputFl: File, signalFl: File) extends CommandLineFunction {
  @Argument(doc = "the bam file name") var bamFile = bamName
  @Input(doc = "the bam to sample file") var bamToSampleFile = bamToSample
  @Output(doc = "the output file") var outputFile = outputFl
  @Input(doc = "the signal file") var signalFile = signalFl
  @Input(doc = "library directory") var libDir = libraryDir
  memoryLimit = Some(4)

  def commandLine = "Rscript %s/R/individual_segmentation.R --bam.file.name %s --bam.to.sample.file %s --output.filename %s --signal.file %s --r.dir %s".format(libDir.getAbsolutePath(), bamFile, bamToSampleFile.getAbsolutePath(), outputFile.getAbsolutePath(), signalFile.getAbsolutePath(),libDir.getAbsolutePath())
}


// post process the data using the R script
class PostProcessData(libraryDir: File, normal_bait_coverage: File, tumor_bait_coverage: File, target_file: File,
                      useCachedData: Boolean, cachedLocation: File, sToRGFileNormals: File, sToRGFileTumors: File,
                      outputDir: File, tangentLocation: File, tangentOutputLocation: File, buildType: String, analysisSet: String,
                      byLane: Boolean, normalBaitFile: File, signalFLs: List[File], signalTSV: File, histoData: Boolean, acovFiles: List[File], sexFl: File) extends CommandLineFunction {
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
  @Input(doc = "information about the sex of every sample") var sexInfo = sexFl

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
  def commandLine = "Rscript %s/R/tangent_normalize.R --normal.lane.data %s --tumor.lane.data %s --target.list %s --use.cache %s --cache.location %s --script.dir %s --normal.sample.to.lanes.file %s --tumor.sample.to.lanes.file %s --output.location %s --tangent.database.location %s --output.tangent.database %s --build %s --analysis.set.name %s --bylane %s --bait.factor %s --signal.files %s --histo.data %s --sex.chromosomes %s".format(libDir.getAbsolutePath(), normalFile.getAbsolutePath(), tumorFile.getAbsolutePath(), targetFile.getAbsolutePath(), useCache, cacheLocation.getAbsolutePath(), libDir.getAbsolutePath(), normalSampleToReadGroupFile.getAbsolutePath(), tumorSampleToReadGroupFile.getAbsolutePath(), outputDirectory.getAbsolutePath(), tangentDatabase.getAbsolutePath(), tangentOutputDatabase.getAbsolutePath(), build, analysisSetName,byLaneData, nbf.getAbsolutePath(), signalTSVFile.getAbsolutePath(), uhd, sexFl.getAbsolutePath())
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

// given a VCF file, split out a BED file for each of the samples into the specified directory
class VCFToBED(vcfFile: File, outputDir: File, normalBamToSample: File, tumorBamToSample: File, individualMapping: File, utilityLoc: File, outputFile: File, sexCalls: File) extends CommandLineFunction {
  @Input(doc = "the VCF input file") var vcf = vcfFile
  @Input(doc = "the output directory to put the BED files in") var out = outputDir
  @Input(doc = "the mapping of normal BAM files to their sample name in BAM file") var nmap = normalBamToSample
  @Input(doc = "the mapping of individual to its normal and tumor bam files") var imap = individualMapping
  @Input(doc = "the mapping of tumor BAM files to their sample names ") var tmap = tumorBamToSample
  @Output(doc = "the final csv output file") var outFile = outputFile
  @Output(doc = "the sex calls for each tumor sample name") var sCalls = sexCalls
  @Argument(doc = "where to find the scripts (we add the utils/python/ to this path)") var loc = utilityLoc
  memoryLimit = Some(2)

  def commandLine = "python %s/utils/python/vcf_to_bed.py --vcf %s --outputdir %s -nmap %s -tmap %s -imap %s -sf %s".format(loc.getAbsolutePath(),vcf.getAbsolutePath(),outputDir.getAbsolutePath(),nmap.getAbsolutePath(),tmap.getAbsolutePath(),imap.getAbsolutePath(), sCalls.getAbsolutePath());
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
class TumorNormalFile(inputFile: String) extends ReadDelimitedFile(inputFile, "\t", true, "sample\ttumor_bam\tnormal_bam\n") {
  var tBamFiles = Map.empty[String, File]
  var nBamFiles = Map.empty[String, File]
  for (line <- rows) {
    var sp = line.split(delim)
    //println(sp.length)
    if (sp.length >= 3 && sp(0) != null && sp(1) != null && sp(2) != null) {
      if (!sp(1).equals("NA") && !sp(1).equals("")) {
        tBamFiles += sp(0) -> new File(sp(1))
        //println("tumor")
      }
      if (!sp(2).equals("NA") && !sp(2).equals("")) {
        nBamFiles += sp(0) -> new File(sp(2))
        //println("normal")
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
}

