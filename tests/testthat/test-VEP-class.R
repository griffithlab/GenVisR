# get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")

test_that("Samples are added from file name", {
    
    # single file
    testFile <- Sys.glob(paste0(testFileDir, "/FLX0070Naive.vep"))
    sample <- as.character(unique(VEP(testFile)@vepObject@sample))
    expect_equal(sample, "FLX0070Naive")
    
    # multiple files
    testFile <- Sys.glob(paste0(testFileDir, "/*vep"))[1:2]
    sample <- nrow(unique(VEP(testFile)@vepObject@sample))
    expect_equal(sample, 2)
})

test_that("Extra columns are properly split", {
    testFile <- Sys.glob(paste0(testFileDir, "/FLX0070Naive.vep"))
    metaFields <- VEP(testFile)@vepObject@meta
    expect_equal(ncol(metaFields), 57)
    
    expect_true("HGNC" %in% metaFields$SYMBOL_SOURCE)
})