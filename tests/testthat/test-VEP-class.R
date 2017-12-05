# get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")
testFile <- Sys.glob(paste0(testFileDir, "/*.vep"))

# define the object for testing
vepObject <- VEP(testFile)

################################################################################
################### test VEP class construction ################################
################################################################################

context("VEP")

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

################################################################################
####################### test accessor methods ##################################
################################################################################

test_that("accessor method getPosition extracts the proper columns", {
    
    expectedCol <- c("Location")
    extractedCol <- colnames(getPosition(vepObject))
    expect_true(all(extractedCol %in% expectedCol))
})

test_that("accessor method getMutation extracts the proper columns", {
    
    expectedCol <- c("Allele", "Consequence")
    extractedCol <- colnames(getMutation(vepObject))
    expect_true(all(extractedCol %in% expectedCol))
})

test_that("accessor method getSample extracts the proper columns", {
    
    expectedCol <- c("sample")
    extractedCol <- colnames(getSample(vepObject))
    expect_true(all(extractedCol %in% expectedCol))
})

test_that("accessor method getMeta extracts all meta columns", {
    
    expectedColNum <- 58
    extractedColNum <- ncol(getMeta(vepObject))
    expect_true(extractedColNum == expectedColNum)
})

################################################################################
############# test the setMutationHierarchy method in Waterfall ################
################################################################################

setMutationHierarchy.out <- suppressWarnings(setMutationHierarchy(vepObject, mutationHierarchy=NULL, verbose=F))
test_that("setMutationHierarchy outputs a data.table with proper columns", {
    
    # test that it is a data.table
    expect_is(setMutationHierarchy.out, "data.table")
    
    # test that it has the proper columns
    actualCol <- colnames(setMutationHierarchy.out)
    expectedCol <- c("mutation", "color", "label")
    expect_true(all(actualCol %in% expectedCol))
})

# define an empty table of mutation hierarchies
emptyMutationHierarchy <- data.table::data.table()

test_that("setMutationHierarchy adds values for missing mutations not specified but in the primary data", {

    # test that a warning message is created
    expect_warning(setMutationHierarchy(vepObject, mutationHierarchy=emptyMutationHierarchy, verbose=FALSE))
    
    # test that output is created for every mutation
    setMutationHierarchy.out <- suppressWarnings(setMutationHierarchy(vepObject, mutationHierarchy=emptyMutationHierarchy, verbose=FALSE))
    expectedMutations <- unique(getMutation(vepObject)$Consequence)
    actualMutations <- setMutationHierarchy.out$mutation
    expect_true(all(expectedMutations %in% actualMutations))
})

# define table with duplicate mutations
duplicateMutationHierarchy <- data.table::data.table("mutation"=c("missense_variant", "missense_variant"), "color"=c("blue", "red"))

test_that("setMutationHierarchy checks for duplicate mutations supplied to input", {
    
    # test that warning is created
    expect_warning(setMutationHierarchy(vepObject, mutationHierarchy=duplicateMutationHierarchy, verbose=FALSE))
    
    # test that the duplicate is removed
    output <- suppressWarnings(setMutationHierarchy(vepObject, mutationHierarchy=duplicateMutationHierarchy, verbose=FALSE)$mutation)
    
    boolean <- !any(duplicated(output))
    expect_true(boolean)
})
