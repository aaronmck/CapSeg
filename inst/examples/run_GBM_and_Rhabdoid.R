# quick comaprison of bait factors, coverage, etc between GBM and Rhabdoid
# 
# Author: aaron
###############################################################################
gbm.base = "/xchip/cga2/aaron/copy_number/svn_R/WES_Segmentation/test_data/GBM/"
rhab.base = "/xchip/cga2/aaron/copy_number/svn_R/WES_Segmentation/test_data/Rhabdoid/"

# GBM data
gbm.normal.cov <- paste(gbm.base,"final.normal.csv",sep="/")
gbm.tumor.cov <- paste(gbm.base,"final.tumor.csv",sep="/")
gbm.normal.samples <- paste(gbm.base,"normalSampleInformation.txt",sep="/")
gbm.tumor.samples <- paste(gbm.base,"tumorSampleInformation.txt",sep="/")
gbm.targets <- 
# Rhabdoid data
rhab.normal.cov <- paste(rhab.base,"final.normal.csv",sep="/")
rhab.tumor.cov <- paste(rhab.base,"final.tumor.csv",sep="/")
rhab.normal.samples <- paste(rhab.base,"normalSampleInformation.txt",sep="/")
rhab.tumor.samples <- paste(rhab.base,"tumorSampleInformation.txt",sep="/")

# source the right files
source("R/command_line_utils.R")
source("R/extract_bait_data.R")
source("R/tangent_normalize.R")

target_list_file = "/xchip/cga2/aaron/copy_number/data/targets/hg18/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly18.baits.csv"
baits <- read.csv(target_list_file,header=TRUE,colClasses=c("character","character","integer","integer"))
baits <- baits[!duplicated(baits$name),]
g.baits <- read.csv("./inst/tests/data/gbm.full.baits.csv",header=TRUE,colClasses=c("character","character","integer","integer"))

processed.data.file <- load.exome.data(gbm.normal.cov,baits,removeBadBaitsAndLanes=F)
cn.hat.data <- cn.hat(processed.data.file)
bait.factor.data <- bait.factor(processed.data.file)
bait.cr <- processed.data.file / (bait.factor.data)