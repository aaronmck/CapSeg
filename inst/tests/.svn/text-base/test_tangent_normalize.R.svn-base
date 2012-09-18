library(testthat)
source("../../R/tangent_normalize.R")
source("../../R/extract_bait_data.R")


# the data sample table to read in
data.file <- "data/samples.table.csv"
baits <- read.csv("data/baits.csv",header=TRUE,colClasses=c("character","character","integer","integer"))
tumorSampleTable = read.table("data/tumorSampleInformation.txt",header=T,stringsAsFactors=F,sep="\t")
normalSampleTable = read.table("data/normalSampleInformation.txt",header=T,stringsAsFactors=F,sep="\t")


# ---------------------------------------------------------------
context("check sample merging")
# ---------------------------------------------------------------
GBM.100.normals.lanes <- load.exome.data("data/GBM.100.baits.normal.lanes.csv",GBM.baits,removeBadBaitsAndLanes=F)
merge.lanes.to.samples(GBM.100.normals.lanes,normalSampleTable)

# ---------------------------------------------------------------
context("tangent normalization process tests")
# ---------------------------------------------------------------

processed.data.file <- load.exome.data(data.file,baits,removeBadBaitsAndLanes=F)
cn.hat.data <- cn.hat(processed.data.file)
bait.factor.data <- bait.factor(processed.data.file)
bait.cr <- processed.data.file / (bait.factor.data)

# load some read data too
GBM.baits <- read.csv("data/gbm.baits.csv",header=TRUE,colClasses=c("character","character","integer","integer"))
GBM.normals <- load.exome.data("data/GBM.200.baits.normal.csv",GBM.baits,removeBadBaitsAndLanes=F)
GBM.tumors <- load.exome.data("data/GBM.200.baits.tumor.csv",GBM.baits,removeBadBaitsAndLanes=F)
exome.data <- tangent.normalize(GBM.tumors,GBM.normals,tumorSampleTable,normalSampleTable,FALSE,print.summary.info=F,optimize.bf=F)

test_that("check that the cn hat data comes out right", {
			# our sums across a sample now should always be 1
			expect_that(sum(cn.hat.data[,1]), equals(1))
			expect_that(sum(cn.hat.data[,8]), equals(1))
		})

test_that("the bait factor comes out looking right (the median is one per bait)", {
			# our medians across a bait should always be 1
			expect_that(median(bait.cr[1,]), equals(1.0))
			expect_that(median(bait.cr[10,]), equals(1.0))
		})

test_that("the tumor and normal for a small segment of GBM come out looking right", {
			# our medians should always be
			expect_that(exome.data,is_a("list"))
			
			expect_that(exome.data$normals.cr,is_a("matrix"))
			expect_that(ncol(exome.data$normals.cr), equals(196))
			expect_that(nrow(exome.data$normals.cr), equals(199))
			
			expect_that(exome.data$tumors.cr,is_a("matrix"))
			expect_that(ncol(exome.data$tumors.cr), equals(194))
			expect_that(nrow(exome.data$tumors.cr), equals(199))			
			expect_that(round(median(exome.data$tumors.cr[,1]), digits = 2), equals(0.97))
			expect_that(round(median(exome.data$tumors.cr[,2]), digits = 2), equals(1.12))			
		})



# test that by lane works as well
GBM.100.baits <- read.csv("data/gbm.baits.100.csv",header=TRUE,colClasses=c("character","character","integer","integer"))
GBM.100.tumors.lanes <- load.exome.data("data/GBM.100.baits.tumor.lanes.csv",GBM.baits,removeBadBaitsAndLanes=F)
exome.data.lanes <- tangent.normalize(GBM.100.tumors.lanes,GBM.100.normals.lanes,tumorSampleTable,normalSampleTable,TRUE,print.summary.info=T,optimize.bf=F)

test_that("the tumor and normal for a small segment of GBM come out looking right", {
			# our medians should always be
			expect_that(exome.data.lanes,is_a("list"))
			expect_that(exome.data.lanes$normals.cr,is_a("matrix"))
			expect_that(exome.data.lanes$tumors.cr,is_a("matrix"))			
		})