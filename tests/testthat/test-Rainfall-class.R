# packges needed
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome <- BSgenome.Hsapiens.UCSC.hg19

# get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")
testFile <- Sys.glob(paste0(testFileDir, "/*.vep"))

# define the object for testing
vepObject <- VEP(testFile)
toRainfall.out <- toRainfall(vepObject, BSgenome=BSgenome, verbose=FALSE)

################################################################################
###### Test Rainfall class and associated functions in constructor #############
################################################################################

context("Rainfall Constructor")

####################### annoRainfall ###########################################
annoRainfall.out <- annoRainfall(toRainfall.out, verbose=FALSE)

test_that("annoRainfall correctly annotates transitions and transversions", {
    
    actual <- as.character(unlist(annoRainfall.out[annoRainfall.out$refAllele == "A" & annoRainfall.out$variantAllele == "G",][1,"trans_tranv"]))
    expected <- "A->G or T->C (TI)"
    expect_equal(actual, expected)
    
    actual <- as.character(unlist(annoRainfall.out[annoRainfall.out$refAllele == "A" & annoRainfall.out$variantAllele == "C",][1,"trans_tranv"]))
    expected <- "A->C or T->G (TV)"
    expect_equal(actual, expected)
    
    actual <- as.character(unlist(annoRainfall.out[annoRainfall.out$refAllele == "T" & annoRainfall.out$variantAllele == "G",][1,"trans_tranv"]))
    expected <- "A->C or T->G (TV)"
    expect_equal(actual, expected)
    
})

test_that("annoRainfall works in verbose mode", {
    
    expect_message(annoRainfall(toRainfall.out, verbose=TRUE))
})

# ####################### calcMutDist ############################################
# calcMutDist.out <- calcMutDist(annoRainfall.out, verbose=FALSE)
# 
# ####################### chrSubset ##############################################
# chrSubset.out <- chrSubset(calcMutDist.out, chromosomes=c("chr8"), verbose=FALSE)
# 
# ####################### annoGenomeCoord ########################################
# annoGenomeCoord.out <- annoGenomeCoord(chrSubset.out, BSgenome=BSgenome, verbose=FALSE)
# 
# ####################### formatSample ###########################################
# formatSample.out <- formatSample(annoGenomeCoord.out, sample="FLX0070Naive", verbose=FALSE)



