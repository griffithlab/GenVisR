# get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")
testFile <- Sys.glob(paste0(testFileDir, "/brca.maf"))

# define the object for testing
mafObject <- MutationAnnotationFormat(testFile)

################################################################################
################## test MutationAnnotationFormat class construction ############
################################################################################

context("MutationAnnotationFormat")

test_that("MutationAnnotationFormat can construct object from a file path", {
    expect_s4_class(mafObject, "MutationAnnotationFormat")
})

test_that("MutationAnnotationFormat stops if version is not supported", {
    expect_error(MutationAnnotationFormat(testFile, version=0))
})

################################################################################
################## test accessor methods #######################################
################################################################################

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

test_that("accessor method getVersion extracts the version number", {
    
    expected <- 2.4
    actual <- getVersion(mafObject)
    expect_equal(expected, actual)
})

test_that("accessor method getPath extracts the appropriate maf file paths", {
    
    expected <- 1
    actual <- length(getPath(mafObject))
    expect_equal(expected, actual)
})

################################################################################
############# test the setMutationHierarchy method in Waterfall ################
################################################################################

setMutationHierarchy.out <- setMutationHierarchy(mafObject, mutationHierarchy=NULL, verbose=F)
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
    expect_warning(setMutationHierarchy(mafObject, mutationHierarchy=emptyMutationHierarchy, verbose=FALSE), "The following mutations")

    # test that output is created for every mutation
    setMutationHierarchy.out <- suppressWarnings(setMutationHierarchy(mafObject, mutationHierarchy=emptyMutationHierarchy, verbose=FALSE))
    expectedMutations <- unique(getMutation(mafObject)$Variant_Classification)
    actualMutations <- setMutationHierarchy.out$mutation
    expect_true(all(expectedMutations %in% actualMutations))
})

# define table with duplicate mutations
duplicateMutationHierarchy <- data.table::data.table("mutation"=c("RNA", "RNA"), "color"=c("blue", "red"))

test_that("setMutationHierarchy checks for duplicate mutations supplied to input", {

    # test that warning is created
    expect_warning(setMutationHierarchy(mafObject, mutationHierarchy=duplicateMutationHierarchy, verbose=FALSE), "was duplicated")

    # test that the duplicate is removed
    output <- suppressWarnings(setMutationHierarchy(mafObject, mutationHierarchy=duplicateMutationHierarchy, verbose=FALSE)$mutation)
    
    boolean <- !any(duplicated(output))
    expect_true(boolean)
})

test_that("setMutationHierarchy errors if the proper columns are not found in hierarchy", {
    mutations <- setMutationHierarchy.out[,c("mutation", "color")]
    colnames(mutations) <- c("wrong", "columns")
    expect_error(setMutationHierarchy(mafObject, mutationHierarchy=mutations, verbose=FALSE))
})

test_that("setMutationHierarchy works in verbose mode", {
    expect_message(setMutationHierarchy(mafObject, mutationHierarchy=NULL, verbose=TRUE))
})

################################################################################
############# test the toWaterfall method in Waterfall #########################
################################################################################

# define additional objects needed for testing
setMutationHierarchy.out <- setMutationHierarchy(mafObject, mutationHierarchy=NULL, verbose=FALSE)
toWaterfall.out <- toWaterfall(mafObject, hierarchy=setMutationHierarchy.out, labelColumn=NULL, verbose=FALSE)

test_that("toWaterfall outputs the correct columns and data types", {

    # check that the data is of the proper class
    expect_is(toWaterfall.out, "data.table")

    # check for the correct columns
    expectedCol <- c("sample", "gene", "mutation", "label")
    actualCol <- colnames(toWaterfall.out)
    expect_true(all(actualCol %in% expectedCol))
})

test_that("toWaterfall adds a specified label column", {

    toWaterfall.out <- toWaterfall(mafObject, hierarchy=setMutationHierarchy.out, labelColumn="Hugo_Symbol", verbose=FALSE)
    expectedValues <- getMeta(mafObject)$Hugo_Symbol
    expect_true(all(toWaterfall.out$label %in% expectedValues))
})

test_that("toWaterfall removes duplicate mutations", {
    
    # create maf object with a duplicate row
    mafObject@mafObject@position <- getPosition(mafObject)[c(1, 1),]
    mafObject@mafObject@mutation <- getMutation(mafObject)[c(1, 1),]
    mafObject@mafObject@sample <- getSample(mafObject)[c(1, 1),]
    mafObject@mafObject@meta <- getMeta(mafObject)[c(1, 1),]

    # create mock waterfall
    toWaterfall.out <- toWaterfall(mafObject, hierarchy=setMutationHierarchy.out, labelColumn=NULL, verbose=FALSE)
    
    expect_true(nrow(toWaterfall.out) == 1)
})

test_that("toWaterfall works in verbose mode", {
    expect_message(toWaterfall(mafObject, hierarchy=setMutationHierarchy.out, labelColumn=NULL, verbose=TRUE))
})

test_that("toWaterfall checks the label parameter", {
    
    expect_warning(toWaterfall(mafObject, hierarchy=setMutationHierarchy.out, labelColumn=c("Mutation_Status", "Validation_Method"), verbose=FALSE))
    expect_warning(toWaterfall(mafObject, hierarchy=setMutationHierarchy.out, labelColumn=c("Not Here"), verbose=FALSE))
})

################################################################################
############# test the toMutSpectra method in MutSpectra #######################
################################################################################

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

test_that("toMutSpectra works in verbose mode", {
    expect_message(toMutSpectra(mafObject, verbose=TRUE))
})

################################################################################
############# test the toRainfall method in Rainfall ###########################
################################################################################

test_that("toRainfall removes duplicate genomic mutations", {
    
    # create vep object with a duplicate row
    mafObject@mafObject@position <- mafObject@mafObject@position[c(1, 1),]
    mafObject@mafObject@mutation <- mafObject@mafObject@mutation[c(1, 1),]
    mafObject@mafObject@sample <- mafObject@mafObject@sample[c(1, 1),]
    mafObject@mafObject@meta <- mafObject@mafObject@meta[c(1, 1),]
    
    toRainfall.out <- suppressWarnings(toRainfall(mafObject, BSgenome=NULL, verbose=FALSE))
    
    expect <- 1
    actual <- nrow(toRainfall.out)
    expect_equal(expect, actual)
})

test_that("toRainfall removes entries with no mutation", {
    
    # create maf object with a row containing no mutation
    mafObject@mafObject@position <- mafObject@mafObject@position[c(1, 1),]
    mafObject@mafObject@mutation <- mafObject@mafObject@mutation[c(1, 1),]
    mafObject@mafObject@sample <- mafObject@mafObject@sample[c(1, 1),]
    mafObject@mafObject@meta <- mafObject@mafObject@meta[c(1, 1),]
    mafObject@mafObject@mutation[1,"Tumor_Seq_Allele2"] <- "T"
    
    toRainfall.out <- suppressWarnings(toRainfall(mafObject, BSgenome=NULL, verbose=FALSE))
    
    expect <- 1
    actual <- nrow(toRainfall.out)
    expect_equal(expect, actual)
})

test_that("toRainfall works in verbose mode", {
    
    expect_message(suppressWarnings(toRainfall(mafObject, BSgenome=NULL, verbose=TRUE)))
})

################################################################################
############# test the toLolliplot method in Lolliplot #########################
################################################################################

toLolliplot.out <- toLolliplot(mafObject, verbose=FALSE)

test_that("toLolliplot outputs the correct columns", {
    
    # test that it has the proper columns
    actualCol <- colnames(toLolliplot.out)
    expectedCol <- c("sample", "chromosome", "start", "stop", "reference", "variant", "gene", "consequence")
    expect_true(all(actualCol %in% expectedCol))
})

test_that("toLolliplot works in verbose mode", {
    
    expect_message(toLolliplot(mafObject, verbose=TRUE))
})

