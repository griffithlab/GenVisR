# get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")
testFile <- Sys.glob(paste0(testFileDir, "/brca.maf"))

# define the object for testing
mafObject <- MutationAnnotationFormat(testFile)

test_that("Vesion number is correctly extracted", {
    
    expectedVersion <- 2.4
    expect_equal(expectedVersion, mafObject@version)
})

test_that("accessor method getPosition extracts the proper columns", {
    
    expectedCol <- c("Chromosome", "Start_Position", "End_Position", "Strand")
    extractedCol <- colnames(getPosition(mafObject))
    expect_true(all(extractedCol %in% expectedCol))
})

test_that("accessor method getMutation extracts the proper columns", {
    
    expectedCol <- c("Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
    extractedCol <- colnames(getMutation(mafObject))
    expect_true(all(extractedCol %in% expectedCol))
})

test_that("accessor method getSample extracts the proper columns", {
    
    expectedCol <- c("Tumor_Sample_Barcode")
    extractedCol <- colnames(getSample(mafObject))
    expect_true(all(extractedCol %in% expectedCol))
})

test_that("accessor method getMeta extracts all meta columns", {
    
    expectedColNum <- 45
    extractedColNum <- ncol(getMeta(mafObject))
    expect_true(extractedColNum == expectedColNum)
})

# # define additional objects needed for testing
# mockWaterfall <- mockS4object("Waterfall")
# mockWaterfall@MutationHierarchy <- setMutationHierarchy(mafObject, mutationHierarchy=NULL, verbose=FALSE)
# 
# test_that("setMutationHierarchy outputs a data.table with proper columns", {
#     
#     # test that it is a data.table
#     expect_is(mockWaterfall@MutationHierarchy, "data.table")
#     
#     # test that it has the proper columns
#     actualCol <- colnames(mockWaterfall@MutationHierarchy)
#     expectedCol <- c("mutation", "color", "label")
#     expect_true(all(actualCol %in% expectedCol))
# })
# 
# # define an empty table of mutation hierarchies
# emptyMutationHierarchy <- data.table::data.table()
# 
# test_that("setMutationHierarchy adds values for missing mutations not specified but in the primary data", {
#     
#     # test that a warning message is created
#     expect_warning(setMutationHierarchy(mafObject, mutationHierarchy=emptyMutationHierarchy, verbose=FALSE), "The following mutations")
#     
#     # test that output is created for every mutation
#     mockWaterfall@MutationHierarchy <- suppressWarnings(setMutationHierarchy(mafObject, mutationHierarchy=emptyMutationHierarchy, verbose=FALSE))
#     expectedMutations <- unique(mafObject@mafObject@mutation$Variant_Classification)
#     actualMutations <- mockWaterfall@MutationHierarchy$mutation
#     expect_true(all(expectedMutations %in% actualMutations))
# })
# 
# # define table with duplicate mutations
# duplicateMutationHierarchy <- data.table::data.table("mutation"=c("RNA", "RNA"), "color"=c("blue", "red"))
# 
# test_that("setMutationHierarchy checks for duplicate mutations supplied to input", {
#     
#     # test that warning is created
#     expect_warning(setMutationHierarchy(mafObject, mutationHierarchy=duplicateMutationHierarchy, verbose=FALSE), "was duplicated")
#     
#     # test that the duplicate is removed
#     output <- suppressWarnings(setMutationHierarchy(mafObject, mutationHierarchy=duplicateMutationHierarchy, verbose=FALSE)$mutation)
#     boolean <- !any(duplicated(output))
#     expect_true(boolean)
# })
# 
# # define additional objects needed for testing
# mockWaterfall@primaryData <- toWaterfall(mafObject, hierarchy=mockWaterfall, labelColumn=NULL, verbose=FALSE)
# 
# test_that("toWaterfall outputs the correct columns and data types", {
#     
#     # check that the data is of the proper class
#     expect_is(mockWaterfall@primaryData, "data.table")
#     
#     # check for the correct columns
#     expectedCol <- c("sample", "gene", "mutation", "label")
#     actualCol <- colnames(mockWaterfall@primaryData)
#     expect_true(all(actualCol %in% expectedCol))
# })
# 
# test_that("toWaterfall adds a specified label column", {
#     
#     mockWaterfall@primaryData <- toWaterfall(mafObject, hierarchy=mockWaterfall, labelColumn="Hugo_Symbol", verbose=FALSE)
#     expectedValues <- mafObject@mafObject@meta$Hugo_Symbol
#     expect_true(all(mockWaterfall@primaryData$label %in% expectedValues))
# })
# 
# test_that("toWaterfall removes duplicate mutations", {
#     
    # # create maf object with a duplicate row
    # mafObject@mafObject@position <- mafObject@mafObject@position[c(1, 1),]
    # mafObject@mafObject@mutation <- mafObject@mafObject@mutation[c(1, 1),]
    # mafObject@mafObject@sample <- mafObject@mafObject@sample[c(1, 1),]
    # mafObject@mafObject@meta <- mafObject@mafObject@meta[c(1, 1),]
#     
#     # create mock waterfall
#     mockWaterfall@primaryData <- toWaterfall(mafObject, hierarchy=mockWaterfall, labelColumn=NULL, verbose=FALSE)
#     
#     expect_true(nrow(mockWaterfall@primaryData) == 1)
# })

# create output to test
primaryData <- toMutSpectra(mafObject, verbose=FALSE)

test_that("toMutSpectra keeps only SNPs", {
    
    boolean <- all(nchar(primaryData$refAllele) == 1 & nchar(primaryData$variantAllele) == 1)
    expect_true(boolean)
})

test_that("toMutSpectra removes duplicate mutations", {
    
    # create maf object with a duplicate row
    mafObject@mafObject@position <- mafObject@mafObject@position[c(1, 1),]
    mafObject@mafObject@mutation <- mafObject@mafObject@mutation[c(1, 1),]
    mafObject@mafObject@sample <- mafObject@mafObject@sample[c(1, 1),]
    mafObject@mafObject@meta <- mafObject@mafObject@meta[c(1, 1),]
    
    # create output to test
    primaryData <- toMutSpectra(mafObject, verbose=FALSE)
    
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
