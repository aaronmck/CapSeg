library(testthat)
source("../../R/command_line_utils.R")
context("command line tests")

# fake command line args
fake.cmd.args = list()
fake.cmd.args$bool.true = "True"
fake.cmd.args$bool.true.weird.case = "TrUe"
fake.cmd.args$bool.false = "False"
fake.cmd.args$bool.false.weird.case = "faLSe"
fake.cmd.args$non.bool = "notabool"
ret = post.process.commandline(fake.cmd.args)

# test that 
test_that("true processing works", {			
			expect_that(ret$bool.true, equals(TRUE))
			expect_that(ret$bool.true.weird.case, equals(TRUE))
		})

test_that("false processing works", {
			expect_that(post.process.commandline(fake.cmd.args)$bool.false, equals(FALSE))
			expect_that(post.process.commandline(fake.cmd.args)$bool.false.weird.case, equals(FALSE))
		})

test_that("normal string processing works", {
			expect_that(post.process.commandline(fake.cmd.args)$non.bool, equals(fake.cmd.args$non.bool))
		})

