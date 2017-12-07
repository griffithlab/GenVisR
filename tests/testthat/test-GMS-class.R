# get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")
testFile <- Sys.glob(paste0(testFileDir, "/FL.gms"))

# define the objects for testing
gmsObject <- GMS(testFile)

################################################################################
##################### test GMS class construction ##############################
################################################################################

context("GMS")

test_that("GMS can construct object from a file path", {
    expect_s4_class(gmsObject, "GMS")
})

test_that("GMS can construct object from data already loaded in R", {
    testData <- fread(testFile)
    expect_s4_class(GMS(data=testData), "GMS")
})

test_that("GMS warns if conversion to a data.table in required", {
    testData <- read.delim(testFile)
    expect_warning(GMS(data=testData))
    expect_s4_class(suppressWarnings(GMS(data=testData)), "GMS")
})

test_that("GMS fills in the sample from the file name if the sample column is absent", {
    gmsObject <- GMS(path=paste0(testFileDir, "/*.gms_test"))
    expect_true(all(as.character(getSample(gmsObject)$sample) %in% c("H_ML-FLX001-FLX001", "H_ML-FLX003-FLX003")))
})

test_that("GMS errors if no files are found", {
    expect_error(GMS(path=paste0(testFileDir, "/*.not_here")))
})

test_that("GMS stops if version is not supported", {
    expect_error(GMS(testFile, version=0))
})
################################################################################
##################### test GMS accessors #######################################
################################################################################

test_that("accessor method getPosition extracts the proper columns", {
    
    expectedCol <- c("chromosome_name", "start", "stop")
    extractedCol <- colnames(getPosition(gmsObject))
    expect_true(all(extractedCol %in% expectedCol))
})

test_that("accessor method getMutation extracts the proper columns", {
    
    expectedCol <- c("reference", "variant", "trv_type")
    extractedCol <- colnames(getMutation(gmsObject))
    expect_true(all(extractedCol %in% expectedCol))
})

test_that("accessor method getSample extracts the proper columns", {
    
    expectedCol <- c("sample")
    extractedCol <- colnames(getSample(gmsObject))
    expect_true(all(extractedCol %in% expectedCol))
})

test_that("accessor method getMeta extracts all meta columns", {
    
    expectedColNum <- 28
    extractedColNum <- ncol(getMeta(gmsObject))
    expect_true(extractedColNum == expectedColNum)
})

test_that("accessor method getVersion extracts the version number", {
    
    expected <- 4
    actual <- getVersion(gmsObject)
    expect_equal(expected, actual)
})

test_that("accessor method getPath extracts the appropriate gms file paths", {
    
    expected <- 1
    actual <- length(getPath(gmsObject))
    expect_equal(expected, actual)
})

################################################################################
############## test the setMutationHierarchy method in Waterfall ###############
################################################################################

setMutationHierarchy.out <- setMutationHierarchy(gmsObject, mutationHierarchy=NULL, verbose=FALSE)
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
    expect_warning(setMutationHierarchy(gmsObject, mutationHierarchy=emptyMutationHierarchy, verbose=FALSE))

    # test that output is created for every mutation
    mutationHierarchy <- suppressWarnings(setMutationHierarchy(gmsObject, mutationHierarchy=emptyMutationHierarchy, verbose=FALSE))
    expectedMutations <- unique(getMutation(gmsObject)$trv_type)
    actualMutations <- mutationHierarchy$mutation
    expect_true(all(expectedMutations %in% actualMutations))
})

# define table with duplicate mutations
duplicateMutationHierarchy <- data.table::data.table("mutation"=c("rna", "rna"), "color"=c("blue", "red"))

test_that("setMutationHierarchy checks for duplicate mutations supplied to input", {
    
    # test that warning is created
    expect_warning(setMutationHierarchy(gmsObject, mutationHierarchy=duplicateMutationHierarchy, verbose=FALSE))
    
    # test that the duplicate is removed
    output <- suppressWarnings(setMutationHierarchy(gmsObject, mutationHierarchy=duplicateMutationHierarchy, verbose=FALSE)$mutation)
    
    boolean <- !any(duplicated(output))
    expect_true(boolean)
})

test_that("setMutationHierarchy warns if input is not a data.table", {
    
    mutations <- as.data.frame(setMutationHierarchy.out[,c("mutation", "color")])
    expect_warning(setMutationHierarchy(gmsObject, mutationHierarchy=mutations, verbose=FALSE))
    
    setMutationHierarchy.out.t3 <- setMutationHierarchy(gmsObject, mutationHierarchy=mutations, verbose=FALSE)
    expect_equivalent(setMutationHierarchy.out.t3, setMutationHierarchy.out)
})

test_that("setMutationHierarchy errors if the proper columns are not found in hierarchy", {
    mutations <- setMutationHierarchy.out[,c("mutation", "color")]
    colnames(mutations) <- c("wrong", "columns")
    expect_error(setMutationHierarchy(gmsObject, mutationHierarchy=mutations, verbose=FALSE))
})

test_that("setMutationHierarchy works in verbose mode", {
    expect_message(setMutationHierarchy(gmsObject, mutationHierarchy=NULL, verbose=TRUE))
})

################################################################################
############# test the toWaterfall method in Waterfall #########################
################################################################################

# define additional objects needed for testing
setMutationHierarchy.out <- setMutationHierarchy(gmsObject, mutationHierarchy=NULL, verbose=FALSE)
toWaterfall.out <- toWaterfall(gmsObject, hierarchy=setMutationHierarchy.out, labelColumn=NULL, verbose=FALSE)

test_that("toWaterfall outputs the correct columns and data types", {
    
    # check that the data is of the proper class
    expect_is(toWaterfall.out, "data.table")
    
    # check for the correct columns
    expectedCol <- c("sample", "gene", "mutation", "label")
    actualCol <- colnames(toWaterfall.out)
    expect_true(all(actualCol %in% expectedCol))
})

test_that("toWaterfall adds a specified label column", {
    toWaterfall.out <- toWaterfall(gmsObject, hierarchy=setMutationHierarchy.out, labelColumn="amino_acid_change", verbose=FALSE)
    expectedValues <- getMeta(gmsObject)$amino_acid_change
    expect_true(all(toWaterfall.out$label %in% expectedValues))
})

test_that("toWaterfall removes duplicate mutations", {
    
    # create maf object with a duplicate row
    gmsObject@gmsObject@position <- getPosition(gmsObject)[c(1, 1),]
    gmsObject@gmsObject@mutation <- getMutation(gmsObject)[c(1, 1),]
    gmsObject@gmsObject@sample <- getSample(gmsObject)[c(1, 1),]
    gmsObject@gmsObject@meta <- getMeta(gmsObject)[c(1, 1),]
    
    # create mock waterfall
    toWaterfall.out <- toWaterfall(gmsObject, hierarchy=setMutationHierarchy.out, labelColumn=NULL, verbose=FALSE)
    
    expect_true(nrow(toWaterfall.out) == 1)
})

test_that("toWaterfall works in verbose mode", {
    expect_message(toWaterfall(gmsObject, hierarchy=setMutationHierarchy.out, labelColumn=NULL, verbose=TRUE))
})

test_that("toWaterfall checks the label parameter", {
    
    expect_warning(toWaterfall(gmsObject, hierarchy=setMutationHierarchy.out, labelColumn=c("VARIANT_CLASS", "BIOTYPE"), verbose=FALSE))
    expect_warning(toWaterfall(gmsObject, hierarchy=setMutationHierarchy.out, labelColumn=c("Not Here"), verbose=FALSE))
})

################################################################################
############# test the toMutSpectra method in MutSpectra #######################
################################################################################

# create output to test
primaryData <- toMutSpectra(gmsObject, verbose=FALSE)

test_that("toMutSpectra keeps only SNPs", {
    
    boolean <- all(nchar(primaryData$refAllele) == 1 & nchar(primaryData$variantAllele) == 1)
    expect_true(boolean)
})

test_that("toMutSpectra removes duplicate mutations", {
    
    # create maf object with a duplicate row
    gmsObject@gmsObject@position <- gmsObject@gmsObject@position[c(1, 1),]
    gmsObject@gmsObject@mutation <- gmsObject@gmsObject@mutation[c(1, 1),]
    gmsObject@gmsObject@sample <- gmsObject@gmsObject@sample[c(1, 1),]
    gmsObject@gmsObject@meta <- gmsObject@gmsObject@meta[c(1, 1),]
    
    # create output to test
    primaryData <- toMutSpectra(gmsObject, verbose=FALSE)
    
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
    expect_message(toMutSpectra(gmsObject, verbose=TRUE))
})
