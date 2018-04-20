# get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")
testFile <- Sys.glob(paste0(testFileDir, "/*.vep"))

library(testthat)
# define the object for testing
path <- "~/OneDrive/Hepatocellular Carcinoma/HCC_SV/Discovery/raw_sv_data/HCC_Manta/*"
data=NULL
version <- "auto"
svCaller <- "Manta"
paired <- "TRUE"
tumorSample <- "H_MU-748892-1209080"
verbose <- FALSE

## Create input for SV 
vcfObject <- extractVariantCallFormat(path=path, data=NULL, version=version, 
                                      svCaller=svCaller, paired=paired, 
                                      tumorColumn=tumorColumn, verbose=verbose)

################################################################################
############# test VariantCallFormat class construction ########################
################################################################################
context("VariantCallFormat")

test_that("VariantCallFormat can construct object from a file path", {
    expect_s4_class(vcfObject, "VariantCallFormat")
})

test_that("VariantCallFormat errors if both ath and data are null", {
    expect_error(extractVariantCallFormat(path=NULL, data=NULL, version="auto", svCaller="fake", paired=TRUE, tumorColumn=11))
})

test_that("VariantCallFormat errors in no files are found", {
    expect_error(extractVariantCallFormat(path=paste0(path, "/*.not_here"), svCaller="Manta"))
})

test_that("VariantCallFormat errors if unsupported sv caller is specified", {
    expect_error(extractVariantCallFormat(path=path, version="auto", svCaller="fake", paired=TRUE, tumorColumn=11))
})

test_that("VariantCallFormat errors if unsupported vcf version is specified", {
    expect_error(extractVariantCallFormat(path=path, version="1.5", svCaller="Manta", paired=TRUE, tumorColumn=11))
})

test_that("VariantCallFormat errors if the data is paired but the paired variable is false", {
    expect_error(extractVariantCallFormat(path=path, version="auto", svCaller="Manta", paired=FALSE, tumorColumn=11))
})

test_that("VariantCallFormat errors if the data is paired but tumorColumn designates a column without sample read support data", {
    expect_error(extractVariantCallFormat(path=path, version="auto", svCaller="Manta", paired=TRUE, tumorColumn=8))
})

test_that("VariantCallFormat can construct object from data already loaded in R", {
    testData <- data.table::fread("~/Google Drive/hcc_sv_dataset.txt")
    expect_s4_class(extractVariantCallFormat(data=testData, version="4.1", svCaller="Manta",
                                             paired=TRUE, tumorColumn=11), "VariantCallFormat")
})

test_that("VariantCallFormat warns if conversion to a data.table is required", {
    testData <- data.table::fread("~/Google Drive/hcc_sv_dataset.txt")
    testData <- as.data.frame(testData)
    expect_warning(extractVariantCallFormat(data=testData, version="4.1", svCaller="Manta",
                                            paired=TRUE, tumorColumn=11))
    expect_s4_class(suppressWarnings(extractVariantCallFormat(data=testData, version="4.1", svCaller="Manta", 
                                                              paired=TRUE, tumorColumn=11)), "VariantCallFormat")
})

test_that("VariantCallFormat errors if data has incorrect columns", {
    testData <- data.table::fread("~/Google Drive/hcc_sv_dataset.txt")
    testData <- testData[,-length(colnames(testData)), with=FALSE]
    expect_error(extractVariantCallFormat(data=testData, version="4.1", svCaller="Manta", 
                                          paired=TRUE, tumorColumn=11))

        
})

test_that("VariantCallFormat errors if data has more than 2 columns that are possibly sample read info", {
    testData <- data.table::fread("~/Google Drive/hcc_sv_dataset.txt")
    testData$fakeData <- "fakeData"
    expect_error(extractVariantCallFormat(data=testData, version="4.1", svCaller="Manta", 
                                          paired=TRUE, tumorColumn=11))
})

