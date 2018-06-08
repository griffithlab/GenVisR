# packges needed
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome <- BSgenome.Hsapiens.UCSC.hg19

# get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")
testFile <- Sys.glob(paste0(testFileDir, "/*.vep"))

# define the object for testing
vepObject <- VEP(testFile)
toRainfall.out <- suppressWarnings(toRainfall(vepObject, BSgenome=BSgenome, verbose=FALSE))

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
calcMutDist.out <- calcMutDist(annoRainfall.out, verbose=FALSE)

test_that("calcMutDist distance calculations reset when jumping between chromosomes", {

    expected <- c(NA, 7.067636, 6.287262)
    actual <- calcMutDist.out[calcMutDist.out$chromosome == "chr20",]
    actual <- actual[order(actual$stop),]
    actual <- actual$log10_difference
    expect_equivalent(expected, actual, tolerance=1e-2)
})

test_that("calcMutDist can calculate distances for duplicate genomic locations", {

    annoRainfall.out <- annoRainfall.out[annoRainfall.out$chromosome == "chr20",]
    annoRainfall.out <- rbind(annoRainfall.out, annoRainfall.out[3,])
    actual <- calcMutDist(annoRainfall.out, verbose=FALSE)
    actual <- actual[actual$start == 60921746]$log10_difference[2]
    expected <- 0
    expect_equivalent(actual, expected)
})

test_that("calcMutDist works in verbose mode", {
    expect_message(calcMutDist(annoRainfall.out, verbose=TRUE))
})

# ####################### chrSubset ##############################################
chrSubset.out <- chrSubset(calcMutDist.out, chromosomes=c("chr8"), verbose=FALSE)

test_that("chrSubset removes all chromosomes not specified", {
    
    actual <- length(unique(chrSubset.out$chromosome))
    expected <- 1
    
    expect_equal(expected, actual)
    
    actual <- chrSubset(calcMutDist.out, chromosomes=c("chr8", "chr1"), verbose=FALSE)
    actual <- length(unique(actual$chromosome))
    expected <- 2
    
    expect_equal(expected, actual)
})

test_that("chrSubset warns if a given chromosome is not found", {
    
    expect_warning(chrSubset(calcMutDist.out, chromosomes=c("chr8", "does_not_exist"), verbose=FALSE))
})

test_that("chrSubset errors if a chromosome subset would eliminate all entries", {
    
    expect_error(expect_warning(chrSubset(calcMutDist.out, chromosomes=c("does_not_exist"), verbose=FALSE)))
})

test_that("chrSubset works in verbose mode", {
    
    expect_message(chrSubset(calcMutDist.out, chromosomes=c("chr8"), verbose=TRUE))
})

####################### annoGenomeCoord ########################################
annoGenomeCoord.out <- annoGenomeCoord(chrSubset.out, BSgenome=BSgenome, verbose=FALSE)

test_that("annoGenomeCoord successfully adds chromosome boundaries from the BSgenome object", {
    
    expected <- 1
    actual <- unique(annoGenomeCoord.out[annoGenomeCoord.out$origin=="chrStart",]$start)
    expect_equal(expected, actual)
    
    actual <- unique(annoGenomeCoord.out[annoGenomeCoord.out$origin=="chrStart",]$stop)
    expect_equal(expected, actual)
    
    # based on chr8
    expected <- 146364022
    actual <- unique(annoGenomeCoord.out[annoGenomeCoord.out$origin=="chrStop",]$start)
    expect_equal(expected, actual)
    
    actual <- unique(annoGenomeCoord.out[annoGenomeCoord.out$origin=="chrStop",]$stop)
    expect_equal(expected, actual)
})

test_that("annoGenomeCoord warns if BSgenome is not set", {
    
    
    expect_warning(annoGenomeCoord(chrSubset.out, BSgenome=NULL, verbose=FALSE))
})

test_that("annoGenomeCoord stops if an unexpected input is given to BSgenome", {
    
    expect_error(annoGenomeCoord(chrSubset.out, BSgenome="character", verbose=FALSE))
})

test_that("annoGenomeCoord works in verbose mode", {
    
    expect_message(annoGenomeCoord(chrSubset.out, BSgenome=BSgenome, verbose=TRUE))
})

# ####################### formatSample ###########################################
formatSample.out <- formatSample(annoGenomeCoord.out, sample="FLX0070Naive", verbose=FALSE)

test_that("formatSample correctly subsets samples if the sample parameter is set", {
    
    expected <- "FLX0070Naive"
    actual <- as.character(unique(formatSample.out$sample))
    expect_equal(expected, actual)
})

