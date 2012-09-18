library(testthat)
source("../../R/extract_bait_data.R")
context("extract data test")

# the data sample table to read in
data.file <- "data/samples.table.csv"

baits <- read.csv("data/baits.csv",header=TRUE,colClasses=c("character","character","integer","integer"))

processed.data.file = load.exome.data(data.file,baits,removeBadBaitsAndLanes=F)
processed.data.file.rbb = load.exome.data(data.file,baits,removeBadBaitsAndLanes=T)

test_that("check that the table is all numeric", {
			# should be 20 rows by 8 columns
			expect_that(nrow(processed.data.file), equals(19))
			expect_that(ncol(processed.data.file), equals(8))
			expect_that(rownames(processed.data.file)[1], equals("bait_1_None"))
			expect_that(rownames(processed.data.file)[10], equals("bait_10_OR4F5"))
			expect_that(rownames(processed.data.file)[19], equals("bait_19_OR4F3"))
			expect_that(median(processed.data.file[,7]), equals(1426.0))
			expect_that(median(processed.data.file[1,]), equals(1620.5))
		})

test_that("check that the table (reduced by a column) is correctly sized", {
			# should be 20 rows by 7 columns
			expect_that(nrow(processed.data.file.rbb), equals(19))
			expect_that(ncol(processed.data.file.rbb), equals(7))
			expect_that(rownames(processed.data.file.rbb)[1], equals("bait_1_None"))
			expect_that(rownames(processed.data.file.rbb)[10], equals("bait_10_OR4F5"))
			expect_that(rownames(processed.data.file.rbb)[19], equals("bait_19_OR4F3"))
			expect_that(median(processed.data.file.rbb[,7]), equals(1426.0))
			expect_that(median(processed.data.file.rbb[1,]), equals(1724))
		})


