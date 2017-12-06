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

# define first test case
setMutationHierarchy.out <- setMutationHierarchy(vepObject, mutationHierarchy=NULL, verbose=FALSE)

# define second test case
extraMutation <- data.table::data.table(mutation=c("missense_variant,splice_region_variant"), color=c("black"))
vepMutations <- data.table::data.table(mutation=setMutationHierarchy.out$mutation, color=setMutationHierarchy.out$color)
vepMutations <- data.table::rbindlist(list(extraMutation, vepMutations), use.names=TRUE, fill=TRUE)
setMutationHierarchy.out.t2 <- setMutationHierarchy(vepObject, mutationHierarchy=vepMutations, verbose=FALSE)

test_that("setMutationHierarchy does not add values for comma delimited vep consequences if they are valid entries", {
    
    # no entires should contain a comma
    expect_true(all(!grepl(",", setMutationHierarchy.out$mutation)))
    
})

test_that("setMutationHierarchy keeps specified comma delimited vep consequences if specifically stated", {

    expect_true(all(vepMutations$mutation %in% setMutationHierarchy.out.t2$mutation))
})

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

################################################################################
############# test the toWaterfall method in Waterfall #########################
################################################################################

# define objects needed for testing 
toWaterfall.out <- suppressWarnings(toWaterfall(vepObject, hierarchy=setMutationHierarchy.out, labelColumn=NULL, verbose=FALSE))

test_that("toWaterfall outputs the correct columns and data types", {
    # check that the data is of the proper class
    expect_is(toWaterfall.out, "data.table")
    
    # check for the correct columns
    expectedCol <- c("sample", "gene", "mutation", "label")
    actualCol <- colnames(toWaterfall.out)
    expect_true(all(actualCol %in% expectedCol))
})

test_that("toWaterfall adds a specified label column", {
    
    toWaterfall.out <- suppressWarnings(toWaterfall(vepObject, hierarchy=setMutationHierarchy.out, labelColumn="VARIANT_CLASS", verbose=FALSE))
    expectedValues <- getMeta(vepObject)$VARIANT_CLASS
    expect_true(all(toWaterfall.out$label %in% expectedValues))
})

test_that("toWaterfall splits VEP supported mutations which are comma delimited and are not stated in the hierarchy", {
    
    ## first test 
    
    # find a case where the mutation should be split
    vepMutations <- unique(getMutation(vepObject)$Consequence)
    vepMutations_shouldbeSplit <- vepMutations[grepl(",", vepMutations, fixed=TRUE) & !vepMutations %in% setMutationHierarchy.out.t2$mutation][7]
    vepMutations_shouldbeSplit_testIndex <- which(as.character(getMutation(vepObject)$Consequence) %in% as.character(vepMutations_shouldbeSplit))[1]
    
    # get the gene and sample for this case
    expected_gene <- getMeta(vepObject)[vepMutations_shouldbeSplit_testIndex]$SYMBOL
    expected_sample <- getSample(vepObject)[vepMutations_shouldbeSplit_testIndex]$sample
    
    # run the waterfall method
    toWaterfall.out <- suppressWarnings(toWaterfall(vepObject, hierarchy=setMutationHierarchy.out.t2, labelColumn=NULL, verbose=FALSE))
    
    # test that for the test case gene/sample the mutation meets the expectation
    actual_test_row <- toWaterfall.out[toWaterfall.out$gene == expected_gene & toWaterfall.out$sample == expected_sample,]
    
    expect_true(actual_test_row$mutation == "missense_variant")
    
    ## second test
    toWaterfall.out <- suppressWarnings(toWaterfall(vepObject, hierarchy=setMutationHierarchy.out, labelColumn=NULL, verbose=FALSE))
    expect_true(all(!grepl(",", toWaterfall.out$mutation, fixed=TRUE)))
})

test_that("toWaterfall does not split comma delimited mutations which are in the mutation hierarcy", {
    
    # find a case where the mutation should not be split
    vepMutations <- unique(getMutation(vepObject)$Consequence)
    vepMutations_shouldNotbeSplit <- vepMutations[grepl(",", vepMutations, fixed=TRUE) & vepMutations %in% setMutationHierarchy.out.t2$mutation][1]
    vepMutations_shouldNotbeSplit_testIndex <- which(as.character(getMutation(vepObject)$Consequence) %in% as.character(vepMutations_shouldNotbeSplit))[1]
    
    # get the gene and sample for this case
    expected_gene <- getMeta(vepObject)[vepMutations_shouldNotbeSplit_testIndex]$SYMBOL
    expected_sample <- getSample(vepObject)[vepMutations_shouldNotbeSplit_testIndex]$sample
    
    # run the waterfall method
    toWaterfall.out <- suppressWarnings(toWaterfall(vepObject, hierarchy=setMutationHierarchy.out.t2, labelColumn=NULL, verbose=FALSE))
    
    # test that for the test case gene/sample the mutation meets the expectation
    actual_test_row <- toWaterfall.out[toWaterfall.out$gene == expected_gene & toWaterfall.out$sample == expected_sample,]
    
    expect_true(actual_test_row$mutation == "missense_variant,splice_region_variant")
})

test_that("toWaterfall properly removes duplicate genomic entries based on a hierarchy", {
    
    # get a test case
    expected_gene <- as.character(getMeta(vepObject)[which(getPosition(vepObject)$Location == "20:60921746"),"SYMBOL"][1])
    expected_sample <- as.character(getSample(vepObject)[which(getPosition(vepObject)$Location == "20:60921746"),"sample"][1])
    
    # run method to test
    toWaterfall.out <- suppressWarnings(toWaterfall(vepObject, hierarchy=setMutationHierarchy.out, labelColumn=NULL, verbose=FALSE))
    
    # test that for the test case gene/sample the mutation meets the expectation
    actual_test_row <- toWaterfall.out[toWaterfall.out$gene == expected_gene & toWaterfall.out$sample == expected_sample,]
    
    expect_equivalent(nrow(actual_test_row), 1)
    expect_true(as.character(actual_test_row$mutation) == "missense_variant")
})

test_that("toWaterfall works in verbose mode", {
    expect_message(suppressWarnings(setMutationHierarchy(vepObject, mutationHierarchy=NULL, verbose=TRUE)))
})

################################################################################
############# test the toMutSpectra method in MutSpectra #######################
################################################################################

library(BSgenome.Hsapiens.UCSC.hg19)

# create output to test
primaryData <- suppressWarnings(toMutSpectra(vepObject, BSgenome=BSgenome.Hsapiens.UCSC.hg19, verbose=FALSE))

test_that("toMutSpectra keeps only SNPs", {
    
    boolean <- all(nchar(primaryData$refAllele) == 1 & nchar(primaryData$variantAllele) == 1)
    expect_true(boolean)
})

test_that("toMutSpectra removes duplicate mutations", {
    
    # create maf object with a duplicate row
    vepObject@vepObject@position <- vepObject@vepObject@position[c(1, 1),]
    vepObject@vepObject@mutation <- vepObject@vepObject@mutation[c(1, 1),]
    vepObject@vepObject@sample <- vepObject@vepObject@sample[c(1, 1),]
    vepObject@vepObject@meta <- vepObject@vepObject@meta[c(1, 1),]
    
    # create output to test
    primaryData <- suppressWarnings(toMutSpectra(vepObject, BSgenome=BSgenome.Hsapiens.UCSC.hg19, verbose=FALSE))
    
    expect_true(nrow(primaryData) == 1)
})

test_that("toMutSpectra creates a data.table with appropriate column names", {
    #test that it is a data.table
    expect_is(primaryData, "data.table")
    
    # test that it has the proper columns
    actualCol <- colnames(primaryData)
    expectedCol <- c("sample", "chromosome", "start", "stop", "refAllele", "variantAllele")
    expect_true(all(actualCol %in% expectedCol))
})

test_that("toMutSpectra grabs the appropriate reference base for a given genomic location", {
    
    # double checked these in ucsc genome browser for hg19
    actual <- primaryData[primaryData$chromosome == "chr1" & primaryData$stop == 2489782,]$refAllele
    expected <- "G"
    expect_equal(actual, expected)
    
    # double checked these in ucsc genome browser for hg19
    actual <- primaryData[primaryData$chromosome == "chr8" & primaryData$stop == 105257280,]$refAllele
    expected <- "A"
    expect_equal(actual, expected)
})

test_that("toMutSpectra works in verbose mode", {
    
    expect_message(toMutSpectra(vepObject, BSgenome=BSgenome.Hsapiens.UCSC.hg19, verbose=TRUE))
})
