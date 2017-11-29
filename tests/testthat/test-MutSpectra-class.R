# get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")
testFile <- Sys.glob(paste0(testFileDir, "/FL.gms"))

# define the objects for testing
gmsObject <- GMS(testFile)
toMutSpectra.out <- toMutSpectra(gmsObject, verbose=FALSE)

################################################################################
## Test MutSpectraPrimaryData class and associated functions in constructor ####
################################################################################

context("MutSpectraPrimaryData Constructor")

##################### test annoMutSpectra ######################################

annoMutSpectra.out <- annoMutSpectra(toMutSpectra.out, verbose=FALSE)

test_that("annoMutSpectra properly classifies reference and variant bases as transitions and transversions", {

    actual <- as.character(unlist(annoMutSpectra.out[annoMutSpectra.out$refAllele == "A" & annoMutSpectra.out$variantAllele == "G",][1,"trans_tranv"]))
    expected <- "A->G or T->C (TI)"
    expect_equal(actual, expected)
    
    actual <- as.character(unlist(annoMutSpectra.out[annoMutSpectra.out$refAllele == "A" & annoMutSpectra.out$variantAllele == "C",][1,"trans_tranv"]))
    expected <- "A->C or T->G (TV)"
    expect_equal(actual, expected)
    
    actual <- as.character(unlist(annoMutSpectra.out[annoMutSpectra.out$refAllele == "T" & annoMutSpectra.out$variantAllele == "G",][1,"trans_tranv"]))
    expected <- "A->C or T->G (TV)"
    expect_equal(actual, expected)
    
})

###################### test calcMutSpectra #####################################

calcMutSpectra.out <- calcMutSpectra(annoMutSpectra.out, verbose=FALSE)

test_that("calcMutSpectra correctly calculates the frequency of transitions/transversions", {

    G2C <- nrow(annoMutSpectra.out[annoMutSpectra.out$refAllele == "G" & annoMutSpectra.out$variantAllele == "C" & annoMutSpectra.out$sample == "FLX007-Naive",])
    C2G <- nrow(annoMutSpectra.out[annoMutSpectra.out$refAllele == "C" & annoMutSpectra.out$variantAllele == "G" & annoMutSpectra.out$sample == "FLX007-Naive",])
    expected <- as.numeric(G2C + C2G)
    actual <- as.numeric(calcMutSpectra.out[calcMutSpectra.out$Sample == "FLX007-Naive" & calcMutSpectra.out$TransTranv == "G->C or C->G (TV)","Frequency"])
    expect_equal(expected, actual)
    
})

test_that("calcMutSpectra correctly calculates the proportion of transitions/transversions", {
    
    G2C <- nrow(annoMutSpectra.out[annoMutSpectra.out$refAllele == "G" & annoMutSpectra.out$variantAllele == "C" & annoMutSpectra.out$sample == "FLX007-Naive",])
    C2G <- nrow(annoMutSpectra.out[annoMutSpectra.out$refAllele == "C" & annoMutSpectra.out$variantAllele == "G" & annoMutSpectra.out$sample == "FLX007-Naive",])
    numerator <- G2C + C2G
    denominator <- nrow(annoMutSpectra.out[annoMutSpectra.out$sample == "FLX007-Naive",])
    expected <- numerator/denominator
    
    actual <- as.numeric(calcMutSpectra.out[calcMutSpectra.out$Sample == "FLX007-Naive" & calcMutSpectra.out$TransTranv == "G->C or C->G (TV)","Proportion"])
    
    expect_equal(expected, actual)
    
})

######################## test sortSamples ######################################

sortSamples.out <- sortSamples(calcMutSpectra.out, sorting="sample", verbose=FALSE)

test_that("sortSamples correctly sorts samples by name", {
    
    expected <- c("FLX001-Naive", "FLX003-Naive", "FLX004-Naive", "FLX005-Naive", "FLX007-Naive")
    actual <- levels(sortSamples.out$Sample)
    
    expect_equal(expected, actual)
    
})

test_that("sortSamples correctly sorts by the observed proportion of mutations", {
    
    sortSamples.out <- sortSamples(calcMutSpectra.out, sorting="mutation", verbose=FALSE)
    expected <- c("FLX007-Naive", "FLX005-Naive", "FLX001-Naive", "FLX004-Naive", "FLX003-Naive")
    actual <- levels(sortSamples.out$Sample)
    
    expect_equal(expected, actual)
})

test_that("sortSamples correctly sorts in a custom order if supplied", {
    
    sampleOrder <- c("FLX005-Naive", "FLX007-Naive", "FLX001-Naive", "FLX004-Naive", "FLX003-Naive")
    sortSamples.out <- sortSamples(calcMutSpectra.out, sorting=sampleOrder, verbose=FALSE)
    expected <- sampleOrder
    actual <- levels(sortSamples.out$Sample)
    
    expect_equal(expected, actual)
    
})

test_that("sortSamples correctly adds samples if not specified in a custom order", {
    
    sampleOrder <- c("FLX007-Naive", "FLX001-Naive", "FLX004-Naive", "FLX003-Naive")
    sortSamples.out <- suppressWarnings(sortSamples(calcMutSpectra.out, sorting=sampleOrder, verbose=FALSE))
    expected <- c(sampleOrder, "FLX005-Naive")
    actual <- levels(sortSamples.out$Sample)
    
    expect_equal(expected, actual)
    
})

test_that("sortSamples correctly removes samples if specified in a custom order", {
    
    sampleOrder <- c("FLX007-Naive", "FLX001-Naive", "FLX004-Naive", "FLX003-Naive", "FLX005-Naive", "not_expected")
    sortSamples.out <- suppressWarnings(sortSamples(calcMutSpectra.out, sorting=sampleOrder, verbose=FALSE))
    expected <- c("FLX007-Naive", "FLX001-Naive", "FLX004-Naive", "FLX003-Naive", "FLX005-Naive")
    actual <- levels(sortSamples.out$Sample)
    
    expect_equal(expected, actual)
    
})
