# packges needed
library(ggplot2)

# get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")
testFile <- Sys.glob(paste0(testFileDir, "/brca.maf"))

# define the objects for testing
mafObject <- MutationAnnotationFormat(testFile)
setMutationHierarchy.out <- setMutationHierarchy(mafObject, mutationHierarchy=NULL, verbose=F)
toWaterfall.out <- toWaterfall(mafObject, hierarchy=setMutationHierarchy.out, labelColumn=NULL, verbose=FALSE)

################################################################################
###### Test WaterfallData class and associated functions in constructor ########
################################################################################

context("WaterfalData Constructor")

#################### setMutationHierarchy ######################################
datatableObject <- data.table::data.table("sample"=rep("test1", 5), "gene"=c(rep("egfr", 2), rep("rb1", 3)),
                                          "mutation"=c(rep("missense", 2), rep("nonsense", 2), "frameshift"),
                                          "amino_acid_change"=rep("p.L546Q", 5))
setMutationHierarchy.out.dt <- suppressWarnings(setMutationHierarchy(datatableObject, mutationHierarchy=NULL, verbose=FALSE))

test_that("setMutationHierarchy outputs a data.table with proper columns", {

    # test that it is a data.table
    expect_is(setMutationHierarchy.out.dt, "data.table")

    # test that it has the proper columns
    actualCol <- colnames(setMutationHierarchy.out.dt)
    expectedCol <- c("mutation", "color", "label")
    expect_true(all(actualCol %in% expectedCol))
})

# define an empty table of mutation hierarchies
emptyMutationHierarchy <- data.table::data.table()

test_that("setMutationHierarchy adds values for missing mutations not specified but in the primary data", {

    # test that a warning message is created
    expect_warning(setMutationHierarchy(datatableObject, mutationHierarchy=emptyMutationHierarchy, verbose=FALSE))

    # test that output is created for every mutation
    mutationHierarchy <- suppressWarnings(setMutationHierarchy(datatableObject, mutationHierarchy=emptyMutationHierarchy, verbose=FALSE))
    expectedMutations <- unique(datatableObject$mutation)
    actualMutations <- mutationHierarchy$mutation
    expect_true(all(expectedMutations %in% actualMutations))
})

# define table with duplicate mutations
duplicateMutationHierarchy <- data.table::data.table("mutation"=c("frameshift", "frameshift"), "color"=c("blue", "red"))

test_that("setMutationHierarchy checks for duplicate mutations supplied to input", {

    # test that warning is created
    expect_warning(setMutationHierarchy(datatableObject, mutationHierarchy=duplicateMutationHierarchy, verbose=FALSE))

    # test that the duplicate is removed
    output <- suppressWarnings(setMutationHierarchy(datatableObject, mutationHierarchy=duplicateMutationHierarchy, verbose=FALSE)$mutation)

    boolean <- !any(duplicated(output))
    expect_true(boolean)
})

test_that("setMutationHierarchy errors if the proper columns are not found in hierarchy", {
    mutations <- setMutationHierarchy.out.dt[,c("mutation", "color")]
    colnames(mutations) <- c("wrong", "columns")
    expect_error(setMutationHierarchy(datatableObject, mutationHierarchy=mutations, verbose=FALSE))
})

test_that("setMutationHierarchy works in verbose mode", {
    expect_message(suppressWarnings(setMutationHierarchy(datatableObject, mutationHierarchy=NULL, verbose=TRUE)))
})

#################### toWaterfall ###############################################

# define additional objects needed for testing
setMutationHierarchy.out.dt <- suppressWarnings(setMutationHierarchy(datatableObject, mutationHierarchy=NULL, verbose=FALSE))
toWaterfall.out.dt <- toWaterfall(datatableObject, hierarchy=setMutationHierarchy.out.dt, labelColumn=NULL, verbose=FALSE)

test_that("toWaterfall outputs the correct columns and data types", {

    # check that the data is of the proper class
    expect_is(toWaterfall.out.dt, "data.table")

    # check for the correct columns
    expectedCol <- c("sample", "gene", "mutation", "label")
    actualCol <- colnames(toWaterfall.out.dt)
    expect_true(all(actualCol %in% expectedCol))

    # also test data frame
    dataframeObject <- as.data.frame(datatableObject)
    toWaterfall.out.df <- toWaterfall(dataframeObject, hierarchy=setMutationHierarchy.out.dt, labelColumn=NULL, verbose=FALSE)
    actualCol <- colnames(toWaterfall.out.df)
    expect_true(all(actualCol %in% expectedCol))
})

test_that("toWaterfall adds a specified label column", {
    toWaterfall.out <- toWaterfall(datatableObject, hierarchy=setMutationHierarchy.out.dt, labelColumn="amino_acid_change", verbose=FALSE)
    expectedValues <- datatableObject$amino_acid_change
    expect_true(all(toWaterfall.out$label %in% expectedValues))
})

test_that("toWaterfall errors if a column is missing", {

    datatableObject <- datatableObject[,c("sample", "gene")]
    expect_error(toWaterfall(datatableObject, hierarchy=setMutationHierarchy.out.dt,  labelColumn=NULL, verbose=FALSE))
})

test_that("toWaterfall checks the label parameter", {

    expect_warning(toWaterfall(datatableObject, hierarchy=setMutationHierarchy.out.dt, labelColumn=c("amino_acid_change", "extra"), verbose=FALSE))
    expect_warning(toWaterfall(datatableObject, hierarchy=setMutationHierarchy.out.dt, labelColumn=c("Not Here"), verbose=FALSE))
})

test_that("toWaterfall works in verbose mode", {

    expect_message(toWaterfall(datatableObject, hierarchy=setMutationHierarchy.out.dt,  labelColumn=NULL, verbose=TRUE))
})

##################### sampSubset ###############################################

test_that("sampSubset filters samples not specified to be kept", {

    testSample <- as.character(toWaterfall.out$sample[1])
    sampSubset.out <- sampSubset(toWaterfall.out, sample=testSample, verbose=F)
    expected <- 1
    actual <- length(unique(sampSubset.out$sample))

    expect_equal(expected, actual)
})

test_that("sampSubset adds samples in if specified and not in the data", {

    testSamples <- c(unique(as.character(toWaterfall.out$sample)), "new_sample")

    expect_warning(sampSubset(toWaterfall.out, sample=testSamples, verbose=F), "but not found")

    sampSubset.out <- suppressWarnings(sampSubset(toWaterfall.out, sample=testSamples, verbose=F))

    expected <- length(unique(toWaterfall.out$sample)) + 1
    actual <- length(unique(sampSubset.out$sample))

    expect_equal(expected, actual)
})

test_that("sampSubset outputs samples as factors", {

    testSample <- as.character(toWaterfall.out$sample[1])
    sampSubset.out <- sampSubset(toWaterfall.out, sample=testSample, verbose=F)

    expect_is(sampSubset.out$sample, "factor")
})

test_that("sampSubset works with verbose mode", {

    testSample <- as.character(toWaterfall.out$sample[1])
    expect_message(sampSubset(toWaterfall.out, samples=testSample, verbose=TRUE))
})

test_that("sampSubset outputs warning if samples is not a character", {

    testSample <- as.factor(toWaterfall.out$sample[1])
    expect_warning(sampSubset(toWaterfall.out, samples=testSample, verbose=FALSE))
})

########################### calcSimpleMutationBurden ###########################

calcSimpleMutationBurden.out <- calcSimpleMutationBurden(toWaterfall.out, coverage=NULL, verbose=FALSE)

test_that("calcSimpleMutationBurden returns only frequencies if coverage is not numeric", {
    expect_true(all(is.na(calcSimpleMutationBurden.out$mutationBurden)))
})

test_that("calcSimpleMutationBurden correctly calculates the frequency of mutations", {
    expected <- nrow(toWaterfall.out[toWaterfall.out$sample == "TCGA-A1-A0SB-01A-11D-A142-09",])
    actual <- as.numeric(calcSimpleMutationBurden.out[calcSimpleMutationBurden.out$sample=="TCGA-A1-A0SB-01A-11D-A142-09", "Freq"])
    expect_equal(expected, actual)
})

test_that("calcSimpleMutationBurden sets samples with no mutation to have a mutation frequency of 0", {

    # check correct number
    testSamples <- c(unique(as.character(toWaterfall.out$sample)), "new_sample")
    sampSubset.out <- suppressWarnings(sampSubset(toWaterfall.out, sample=testSamples, verbose=F))
    calcSimpleMutationBurden.out <- calcSimpleMutationBurden(sampSubset.out, coverage=NULL, verbose=FALSE)
    expected <- 0
    actual <- as.numeric(calcSimpleMutationBurden.out[calcSimpleMutationBurden.out$sample=="new_sample", "Freq"])
    expect_equal(expected, actual)

    # check no duplicate rows
    expected <- 1
    actual <- nrow(calcSimpleMutationBurden.out[calcSimpleMutationBurden.out$sample=="new_sample", "Freq"])
    expect_equal(expected, actual)
})

test_that("calcSimpleMutationBurden calculates a mutation burden if coverage is given", {
    calcSimpleMutationBurden.out <- calcSimpleMutationBurden(toWaterfall.out, coverage=1000000, verbose=FALSE)
    expected <- nrow(toWaterfall.out[toWaterfall.out$sample == "TCGA-A1-A0SB-01A-11D-A142-09",])
    actual <- as.numeric(calcSimpleMutationBurden.out[calcSimpleMutationBurden.out$sample=="TCGA-A1-A0SB-01A-11D-A142-09", "mutationBurden"])
    expect_equal(expected, actual)
})

test_that("calcSimpleMutationBurden works in verbose mode", {

    expect_message(calcSimpleMutationBurden(toWaterfall.out, coverage=NULL, verbose=TRUE))
})

test_that("calcSimpleMutationBurden warns if coverage does not have values for all samples", {

    coverageValues <- c("TCGA-A1-A0SF-01A-11D-A142-09"=5, "TCGA-A1-A0SD-01A-11D-A10Y-09"=4,
                        "TCGA-A1-A0SG-01A-11D-A142-09"=2, "TCGA-A1-A0SE-01A-11D-A099-09"=3)
    expect_warning(calcSimpleMutationBurden(toWaterfall.out, coverage=coverageValues, verbose=FALSE))
})

test_that("calcSimpleMutationBurden warns if coverage contains duplicate names", {
    coverageValues <- c("TCGA-A1-A0SF-01A-11D-A142-09"=5, "TCGA-A1-A0SD-01A-11D-A10Y-09"=4,
                        "TCGA-A1-A0SG-01A-11D-A142-09"=2, "TCGA-A1-A0SE-01A-11D-A099-09"=3,
                        "TCGA-A1-A0SB-01A-11D-A142-09"=1, "TCGA-A1-A0SB-01A-11D-A142-09"=2)
    expect_warning(calcSimpleMutationBurden(toWaterfall.out, coverage=coverageValues, verbose=FALSE))
})

test_that("calcSimpleMutationBurden warns if coverage contains an unexpected name", {
    coverageValues <- c("TCGA-A1-A0SF-01A-11D-A142-09"=5, "TCGA-A1-A0SD-01A-11D-A10Y-09"=4,
                        "TCGA-A1-A0SG-01A-11D-A142-09"=2, "TCGA-A1-A0SE-01A-11D-A099-09"=3,
                        "TCGA-A1-A0SB-01A-11D-A142-09"=1, "extra"=1)
    expect_warning(calcSimpleMutationBurden(toWaterfall.out, coverage=coverageValues, verbose=FALSE))
})

test_that("calcSimpleMutationBurden warns if coverage is not numeric", {

    testCov <- as.character(10000)
    expect_warning(calcSimpleMutationBurden(toWaterfall.out, coverage=testCov, verbose=FALSE))
})

############################ calcComplexMutationBurden #########################

calcComplexMutationBurden.out <- calcComplexMutationBurden(toWaterfall.out, coverage=NULL, verbose=FALSE)

test_that("calcComplexMutationBurden returns only frequencies if coverage is not numeric", {
    expect_true(all(is.na(calcComplexMutationBurden.out$mutationBurden)))
})

test_that("calcComplexMutationBurden correctly calculates the frequencies of mutations", {

    actual <- as.numeric(calcComplexMutationBurden.out[calcComplexMutationBurden.out$sample=="TCGA-A1-A0SB-01A-11D-A142-09" & calcComplexMutationBurden.out$mutation=="Nonsense_Mutation","Freq"])
    expected <- 0
    expect_equal(actual, expected)

    actual <- as.numeric(calcComplexMutationBurden.out[calcComplexMutationBurden.out$sample=="TCGA-A1-A0SB-01A-11D-A142-09" & calcComplexMutationBurden.out$mutation=="Missense_Mutation","Freq"])
    expected <- 11
    expect_equal(actual, expected)
})

test_that("calcComplexMutationBurden sets samples with no mutation to have a mutation frequency of 0", {

    # check correct number
    testSamples <- c(unique(as.character(toWaterfall.out$sample)), "new_sample")
    sampSubset.out <- suppressWarnings(sampSubset(toWaterfall.out, sample=testSamples, verbose=F))
    calcComplexMutationBurden.out <- calcComplexMutationBurden(sampSubset.out, coverage=NULL, verbose=FALSE)
    expected <- 0
    actual <- as.numeric(calcComplexMutationBurden.out[calcComplexMutationBurden.out$sample=="new_sample" & calcComplexMutationBurden.out$mutation=="Missense_Mutation", "Freq"])
    expect_equal(expected, actual)

    # check no duplicate rows (i.e. one row for each mutation type)
    expected <- 18
    actual <- nrow(calcComplexMutationBurden.out[calcComplexMutationBurden.out$sample=="new_sample", "Freq"])
    expect_equal(expected, actual)
})

test_that("calcSimpleMutationBurden calculates a mutation burden if coverage is given", {
    calcComplexMutationBurden.out <- calcComplexMutationBurden(toWaterfall.out, coverage=1000000, verbose=FALSE)
    expected <- 11
    actual <- as.numeric(calcComplexMutationBurden.out[calcComplexMutationBurden.out$sample=="TCGA-A1-A0SB-01A-11D-A142-09" & calcComplexMutationBurden.out$mutation=="Missense_Mutation", "mutationBurden"])
    expect_equal(expected, actual)
})

test_that("calcComplexMutationBurden works in verbose mode", {
    expect_message(calcComplexMutationBurden(toWaterfall.out, coverage=NULL, verbose=TRUE))
})

test_that("calcComplexMutationBurden warns if coverage does not have values for all samples", {

    coverageValues <- c("TCGA-A1-A0SF-01A-11D-A142-09"=5, "TCGA-A1-A0SD-01A-11D-A10Y-09"=4,
                        "TCGA-A1-A0SG-01A-11D-A142-09"=2, "TCGA-A1-A0SE-01A-11D-A099-09"=3)
    expect_warning(calcComplexMutationBurden(toWaterfall.out, coverage=coverageValues, verbose=FALSE))
})

test_that("calcComplexMutationBurden warns if coverage contains duplicate names", {
    coverageValues <- c("TCGA-A1-A0SF-01A-11D-A142-09"=5, "TCGA-A1-A0SD-01A-11D-A10Y-09"=4,
                        "TCGA-A1-A0SG-01A-11D-A142-09"=2, "TCGA-A1-A0SE-01A-11D-A099-09"=3,
                        "TCGA-A1-A0SB-01A-11D-A142-09"=1, "TCGA-A1-A0SB-01A-11D-A142-09"=2)
    expect_warning(calcComplexMutationBurden(toWaterfall.out, coverage=coverageValues, verbose=FALSE))
})

test_that("calcComplexMutationBurden warns if coverage contains an unexpected name", {
    coverageValues <- c("TCGA-A1-A0SF-01A-11D-A142-09"=5, "TCGA-A1-A0SD-01A-11D-A10Y-09"=4,
                        "TCGA-A1-A0SG-01A-11D-A142-09"=2, "TCGA-A1-A0SE-01A-11D-A099-09"=3,
                        "TCGA-A1-A0SB-01A-11D-A142-09"=1, "extra"=1)
    expect_warning(calcComplexMutationBurden(toWaterfall.out, coverage=coverageValues, verbose=FALSE))
})

test_that("calcComplexMutationBurden warns if coverage is not numeric", {

    testCov <- as.character(10000)
    expect_warning(calcComplexMutationBurden(toWaterfall.out, coverage=testCov, verbose=FALSE))
})

############################### rmvMutation ####################################

test_that("rmvMutation warns if a mutation was specified to be kept but is not found", {
    expect_warning(rmvMutation(toWaterfall.out, mutation=c("invalid"), verbose=FALSE), "were not found")
})

rmvMutation.out <- rmvMutation(toWaterfall.out, mutation=c("Missense_Mutation"), verbose=FALSE)
test_that("rmvMutation keeps only mutation specified to be kept", {
    expect_true(all(rmvMutation.out$mutation == "Missense_Mutation"))
})

test_that("rmvMutation works in verbose mode", {
    expect_message(rmvMutation(toWaterfall.out, mutation="Missense_Mutation", verbose=TRUE))
})

test_that("rmvMutation outputs warning if mutation is not of class character", {

    expect_warning(rmvMutation(toWaterfall.out, mutation=as.factor("Missense_Mutation"), verbose=FALSE))
})

############################## mutHiearchySubset ###############################

mutHierarchySubset.out <- mutHierarchySubset(toWaterfall.out, mutationHierarchy=setMutationHierarchy.out, verbose=FALSE)

test_that("mutHierarchySubset removes duplicates entries for a gene/sample in the proper order", {

    # initial test
    expected <- "Missense_Mutation"
    actual <- as.character(mutHierarchySubset.out[mutHierarchySubset.out$sample=="TCGA-A1-A0SF-01A-11D-A142-09" & mutHierarchySubset.out$gene=="MUC5B",]$mutation)
    expect_equal(expected, actual)

    # secondary test, check reverse
    mutHierarchySubset.out <- mutHierarchySubset(toWaterfall.out, mutationHierarchy=setMutationHierarchy.out[rev(1:nrow(setMutationHierarchy.out)),], verbose=FALSE)
    expected <- "Silent"
    actual <- as.character(mutHierarchySubset.out[mutHierarchySubset.out$sample=="TCGA-A1-A0SF-01A-11D-A142-09" & mutHierarchySubset.out$gene=="MUC5B",]$mutation)
    expect_equal(expected, actual)
})

test_that("mutHierarchySubset leaves on one entry per gene/sample", {
    expect_true(all(table(mutHierarchySubset.out$sample, mutHierarchySubset.out$gene) <= 1))
})

test_that("mutHierarchySubset works in verbose mode", {
    expect_message(mutHierarchySubset(toWaterfall.out, mutationHierarchy=setMutationHierarchy.out, verbose=TRUE))
})


############################## geneSubset ######################################

geneSubset.out <- geneSubset(mutHierarchySubset.out, genes=NULL, verbose=FALSE)
test_that("geneSubset keeps samples by outputing an NA value in addition to genes", {

    expect_true(any(is.na(geneSubset.out)))

    geneSubset.out <- geneSubset(mutHierarchySubset.out, genes=c("ARC"), verbose=FALSE)
    expect_true(any(is.na(geneSubset.out)))
})

test_that("geneSubset works in verbose mode", {
    expect_message(geneSubset(mutHierarchySubset.out, genes=c("ARC"), verbose=TRUE))
})

test_that("geneSubset outputs warning if genes in not of class character", {

    expect_warning(geneSubset(mutHierarchySubset.out, genes=as.factor(c("ARC")), verbose=FALSE))
})

############################## recurrenceSubset ################################

recurrenceSubset.out <- recurrenceSubset(mutHierarchySubset.out, recurrence=.4, verbose=FALSE)
test_that("recurrenceSubset correctly outputs a vector of genes above a set recurrence", {
    geneProp <- table(mutHierarchySubset.out[,c("gene")])/length(unique(mutHierarchySubset.out$sample))
    expected <- c(names(geneProp[geneProp >= .4]), NA)
    actual <- recurrenceSubset.out
    expect_true(all(actual %in% expected))
})

test_that("recurrenceSubset correctly resets the recurrence parameter if above the max observed recurrence", {
    expect_warning(recurrenceSubset(mutHierarchySubset.out, recurrence=.5, verbose=FALSE), "exceeds the recurrence")

    recurrenceSubset.out <- suppressWarnings(recurrenceSubset(mutHierarchySubset.out, recurrence=.8, verbose=FALSE))
    geneProp <- table(mutHierarchySubset.out[,c("gene")])/length(unique(mutHierarchySubset.out$sample))
    expected <- c(names(geneProp[geneProp >= .4]), NA)
    actual <- recurrenceSubset.out
    expect_true(all(actual %in% expected))
})

test_that("recurrenceSubset works in verbose mode", {
    expect_message(recurrenceSubset(mutHierarchySubset.out, recurrence=.4, verbose=TRUE))
})

test_that("recurrenceSubset warns if recurrence is not of class numeric", {

    expect_warning(recurrenceSubset(mutHierarchySubset.out, recurrence=as.character(.4), verbose=FALSE))
})

test_that("recurrenceSubset warns if recurrence contains more than one element", {

    expect_warning(recurrenceSubset(mutHierarchySubset.out, recurrence=c(.4, .2), verbose=FALSE))
})

######################### geneFilter ###########################################

keepGenes <- unique(c(geneSubset.out, recurrenceSubset.out))
geneFilter.out <- geneFilter(mutHierarchySubset.out, genes=keepGenes, verbose=FALSE)

test_that("geneFilter correctly subsets data to only genes specified", {
    expect_true(all(geneFilter.out$gene %in% keepGenes))
})

test_that("geneFilter correctly identifies if a gene is specified to be filtered but is not found", {
    keepGenes <- c(keepGenes, "not_in_data")
    expect_warning(geneFilter(mutHierarchySubset.out, genes=keepGenes, verbose=FALSE), "not found in")

    geneFilter.out <- suppressWarnings(geneFilter(mutHierarchySubset.out, genes=keepGenes, verbose=FALSE))
    expect_true("not_in_data" %in% geneFilter.out$gene)
})


test_that("geneFilter works in verbose mode", {
    keepGenes <- unique(c(geneSubset.out, recurrenceSubset.out))
    expect_message(geneFilter(mutHierarchySubset.out, genes=keepGenes, verbose=TRUE))
})

############################# orderGenes #######################################

orderGenes.out <- orderGenes(geneFilter.out, geneOrder=NULL, verbose=FALSE)
test_that("orderGenes correctly refactors genes based on observed frequencies if no gene order is given", {
    geneFreq <- table(orderGenes.out$gene)
    maxGene <- geneFreq[geneFreq == max(geneFreq)]
    maxGene <- unlist(as.list(maxGene))
    expected <- tail(levels(orderGenes.out$gene), n=1) %in% names(maxGene)
    expect_true(expected)
})

test_that("orderGenes correctly refactors genes based on input to the geneOrder parameter", {
    orderGenes.out <- orderGenes(geneFilter.out, geneOrder=c("FAM182B", "CECR2"), verbose=FALSE)
    expected <- all(tail(levels(orderGenes.out), n=2) == c("FAM182B", "CECR2"))
    expect_true(expected)
})

test_that("orderGenes handles cases where genes are supplied which are not in the data", {
    expect_warning(orderGenes(geneFilter.out, geneOrder=c("FAM182B", "CECR2", "notPresent"), verbose=FALSE), "The following arguments to geneOrder were")

    # check that an order is still set correctly
    orderGenes.out <- suppressWarnings(orderGenes(geneFilter.out, geneOrder=c("FAM182B", "CECR2", "notPresent"), verbose=FALSE))
    expected <- all(tail(levels(orderGenes.out), n=2) == c("FAM182B", "CECR2"))
    expect_true(expected)

    # check that the missing gene is not included in the levels
    expect_true(!"notPresent" %in% levels(orderGenes.out$gene))
})

test_that("orderGenes detects if no genes supplied are present in the data and acts accordingly", {
    expect_warning(orderGenes(geneFilter.out, geneOrder=c("notPresent"), verbose=FALSE), "Found no genes")

    orderGenes.out <- suppressWarnings(orderGenes(geneFilter.out, geneOrder=c("notPresent"), verbose=FALSE))
    geneFreq <- table(orderGenes.out$gene)
    maxGene <- geneFreq[geneFreq == max(geneFreq)]
    maxGene <- unlist(as.list(maxGene))
    expected <- tail(levels(orderGenes.out$gene), n=1) %in% names(maxGene)
    expect_true(expected)
})

test_that("orderGenes works in verbose mode", {
    expect_message(orderGenes(geneFilter.out, geneOrder=NULL, verbose=TRUE))
})

test_that("orderGenes warns if geneOrder is not of class character", {

   expect_warning(orderGenes(geneFilter.out, geneOrder=as.factor(c("FAM182B", "CECR2")), verbose=FALSE))
})

test_that("orderGenes warns if geneOrder contains duplicates", {

    expect_warning(orderGenes(geneFilter.out, geneOrder=c("FAM182B", "CECR2", "FAM182B"), verbose=FALSE))
})

############################# maxGeneSubset ####################################

maxGeneSubset.out <- maxGeneSubset(orderGenes.out, geneMax=2, verbose=FALSE)

test_that("maxGeneSubset correctly limits genes to a maximum number", {

    # test1
    expected <- 2
    actual <- length(unique(maxGeneSubset.out$gene))
    expect_equal(expected, actual)

    #test2
    maxGeneSubset.out <- maxGeneSubset(orderGenes.out, geneMax=1, verbose=FALSE)
    expected <- 1
    actual <- length(unique(maxGeneSubset.out$gene))
    expect_equal(expected, actual)
})

test_that("maxGeneSubset correctly deals with situations where geneMax is not an integer", {

    expect_warning(maxGeneSubset(orderGenes.out, geneMax=1.3, verbose=FALSE), "not a whole number")
    maxGeneSubset.out <- suppressWarnings(maxGeneSubset(orderGenes.out, geneMax=1.3, verbose=FALSE))
    expected <- 1
    actual <- length(unique(maxGeneSubset.out$gene))
    expect_equal(expected, actual)
})

test_that("maxGeneSubset works in verbose mode", {
    expect_message(maxGeneSubset(orderGenes.out, geneMax=2, verbose=TRUE))
})

test_that("maxGeneSubset warns if geneMax is not numeric", {
    expect_warning(maxGeneSubset(orderGenes.out, geneMax=as.character(2), verbose=FALSE))
})

############################# orderSamples #####################################

orderSamples.out <- orderSamples(maxGeneSubset.out, sampleOrder=NULL, verbose=FALSE)

test_that("orderSamples correctly reorders samples based on a hierarchy", {

    # test1
    expected1 <- c("TCGA-A1-A0SF-01A-11D-A142-09", "TCGA-A1-A0SD-01A-11D-A10Y-09",
                   "TCGA-A1-A0SG-01A-11D-A142-09", "TCGA-A1-A0SE-01A-11D-A099-09",
                   "TCGA-A1-A0SB-01A-11D-A142-09")
    expected2 <- c("TCGA-A1-A0SF-01A-11D-A142-09", "TCGA-A1-A0SD-01A-11D-A10Y-09",
                   "TCGA-A1-A0SG-01A-11D-A142-09", "TCGA-A1-A0SB-01A-11D-A142-09",
                   "TCGA-A1-A0SE-01A-11D-A099-09")
    actual <- levels(orderSamples.out$sample)
    expect_true(any(c(all(expected1 == actual), all(expected2 == actual))))

    # test2 special case is needed for this to ensure proper ordering, this is created here
    testDT <- data.table::data.table("sample"=c("testSamp1", "testSamp1", "testSamp2", "testSamp3"),
                                     "gene"=c("FAM182B", "ZNF91", "ZNF91", "ZNF740"),
                                     "mutation"=rep("RNA", 4),
                                     "label"=NA,
                                     "count"=2)
    testDT$sample <- factor(testDT$sample, levels=c(levels(testDT$sample), "testSamp1", "testSamp2", "testSamp3"))
    testDT <- data.table::rbindlist(list(maxGeneSubset.out, testDT))
    orderSamples.out <- orderSamples(testDT, sampleOrder=NULL, verbose=FALSE)
    # both outcomes are correct
    expected1 <- c("TCGA-A1-A0SF-01A-11D-A142-09", "testSamp1", "TCGA-A1-A0SD-01A-11D-A10Y-09",
                   "TCGA-A1-A0SG-01A-11D-A142-09", "testSamp2", "testSamp3",
                   "TCGA-A1-A0SB-01A-11D-A142-09", "TCGA-A1-A0SE-01A-11D-A099-09")
    expected2 <- c("TCGA-A1-A0SF-01A-11D-A142-09", "testSamp1", "TCGA-A1-A0SD-01A-11D-A10Y-09",
                   "TCGA-A1-A0SG-01A-11D-A142-09", "testSamp2", "testSamp3",
                   "TCGA-A1-A0SE-01A-11D-A099-09", "TCGA-A1-A0SB-01A-11D-A142-09")
    actual <- levels(orderSamples.out$sample)
    expect_true(any(c(all(expected1 == actual), all(expected2 == actual))))
})

test_that("orderSamples correctly orders samples in a custom order if specified", {

    sampleCustomOrder <- c("TCGA-A1-A0SF-01A-11D-A142-09", "TCGA-A1-A0SD-01A-11D-A10Y-09",
                           "TCGA-A1-A0SG-01A-11D-A142-09", "TCGA-A1-A0SE-01A-11D-A099-09",
                           "TCGA-A1-A0SB-01A-11D-A142-09")
    orderSamples.out <- orderSamples(maxGeneSubset.out, sampleOrder=sampleCustomOrder, verbose=FALSE)
    expected <- sampleCustomOrder
    actual <- levels(orderSamples.out$sample)
    expect_true(all(expected == actual))
})

test_that("orderSamples removes any duplicate samples provided when using a custom order", {

    sampleCustomOrder <- c("TCGA-A1-A0SF-01A-11D-A142-09", "TCGA-A1-A0SD-01A-11D-A10Y-09",
                           "TCGA-A1-A0SG-01A-11D-A142-09", "TCGA-A1-A0SE-01A-11D-A099-09",
                           "TCGA-A1-A0SB-01A-11D-A142-09", "TCGA-A1-A0SF-01A-11D-A142-09")
    expect_warning(orderSamples(maxGeneSubset.out, sampleOrder=sampleCustomOrder, verbose=FALSE), "Found duplicate elements")

    orderSamples.out <- suppressWarnings(orderSamples(maxGeneSubset.out, sampleOrder=sampleCustomOrder, verbose=FALSE))
    expected <- sampleCustomOrder[1:5]
    actual <- levels(orderSamples.out$sample)
    expect_equal(expected, actual)
})

test_that("orderSamples adds extra samples to the data if provided", {

    # check that warning is produced
    sampleCustomOrder <- c("TCGA-A1-A0SF-01A-11D-A142-09", "TCGA-A1-A0SD-01A-11D-A10Y-09",
                           "TCGA-A1-A0SG-01A-11D-A142-09", "TCGA-A1-A0SE-01A-11D-A099-09",
                           "TCGA-A1-A0SB-01A-11D-A142-09", "new_sample")
    expect_warning(orderSamples(maxGeneSubset.out, sampleOrder=sampleCustomOrder, verbose=FALSE), "were not detected")

    # check sample is added to levels
    orderSamples.out <- suppressWarnings(orderSamples(maxGeneSubset.out, sampleOrder=sampleCustomOrder, verbose=FALSE))
    expected <- sampleCustomOrder
    actual <- levels(orderSamples.out$sample)
    expect_equal(expected, actual)

    # check sample is added to data
    expect_true("new_sample" %in% orderSamples.out$sample)
})

test_that("orderSamples works in verbose mode", {
    expect_message(orderSamples(maxGeneSubset.out, sampleOrder=NULL, verbose=TRUE))
})

test_that("orderSamples warns if sampleOrder is not a character", {
    test_sampleOrder <- as.list(c("TCGA-A1-A0SF-01A-11D-A142-09", "TCGA-A1-A0SD-01A-11D-A10Y-09",
                          "TCGA-A1-A0SG-01A-11D-A142-09", "TCGA-A1-A0SE-01A-11D-A099-09",
                          "TCGA-A1-A0SB-01A-11D-A142-09"))
    expect_warning(orderSamples(maxGeneSubset.out, sampleOrder=test_sampleOrder, verbose=TRUE))
})

test_that("orderSamples warns if sampleOrder contains a sample not in the primary data", {

    test_sampleOrder <- c("TCGA-A1-A0SF-01A-11D-A142-09", "TCGA-A1-A0SD-01A-11D-A10Y-09",
                          "TCGA-A1-A0SG-01A-11D-A142-09", "TCGA-A1-A0SE-01A-11D-A099-09",
                          "TCGA-A1-A0SB-01A-11D-A142-09", "does_not_belong")
    expect_warning(orderSamples(maxGeneSubset.out, sampleOrder=test_sampleOrder, verbose=TRUE))

    actual <- suppressWarnings(orderSamples(maxGeneSubset.out, sampleOrder=test_sampleOrder, verbose=TRUE))
    expect_true(!"does_not_belong" %in% actual)
})

############################# constructGeneData #################################
constructGeneData.out <- constructGeneData(orderSamples.out, verbose=FALSE)

test_that("constructGeneData properly summarizes gene reccurrence", {
    expected <- 1
    actual <- constructGeneData.out[constructGeneData.out$gene == "FAM182B" & constructGeneData.out$mutation == "Missense_Mutation",]$count
    expect_equal(expected, actual)

    expected <- 2
    actual <- constructGeneData.out[constructGeneData.out$gene == "CECR2" & constructGeneData.out$mutation == "Missense_Mutation",]$count
    expect_equal(expected, actual)
})

test_that("constructGeneData works in verbose mode", {
    expect_message(constructGeneData(orderSamples.out, verbose=TRUE))
})

######### test WaterfallData class construction with various parameters ########

WaterfallData.out <- WaterfallData(mafObject, labelColumn=NULL, samples=NULL,
                                   mutationHierarchy=NULL, coverage=44000000, mutation=NULL,
                                   genes=NULL, recurrence=NULL, geneOrder=NULL, geneMax=15,
                                   sampleOrder=NULL, verbose=FALSE)

test_that("WaterfallData constructor outputs a S4 class object", {
    expect_s4_class(WaterfallData.out, "WaterfallData")
})

################################################################################
######### test the WaterfallPlots constructor and various methods ##############
################################################################################

################## buildMutationPlot ###########################################

context("Waterfall Mutation Plot")

test_that("buildMutationPlot draws a complex frequency plot correctly", {

    buildMutationPlot.out <- buildMutationPlot(WaterfallData.out, plotA="frequency", plotATally="complex", plotALayers=NULL, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("mutation plot frequency complex", grid::grid.draw(buildMutationPlot.out))

})

test_that("buildMutationPlot draws a simple frequency plot correctly", {

    buildMutationPlot.out <- buildMutationPlot(WaterfallData.out, plotA="frequency", plotATally="simple", plotALayers=NULL, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("mutation plot frequency simple", grid::grid.draw(buildMutationPlot.out))

})

test_that("buildMutationPlot draws a simple burden plot correctly", {

    buildMutationPlot.out <- buildMutationPlot(WaterfallData.out, plotA="burden", plotATally="simple", plotALayers=NULL, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("mutation plot burden simple", grid::grid.draw(buildMutationPlot.out))

})

test_that("buildMutationPlot draws a complex burden plot correctly", {

    buildMutationPlot.out <- buildMutationPlot(WaterfallData.out, plotA="burden", plotATally="complex", plotALayers=NULL, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("mutation plot burden complex", grid::grid.draw(buildMutationPlot.out))

})

test_that("buildMutationPlot works in verbose mode", {

    expect_message(buildMutationPlot(WaterfallData.out, plotA="burden", plotATally="complex", plotALayers=NULL, verbose=TRUE))
})

test_that("buildMutationPlot warns if input to plotA is unexpected", {

    expect_warning(buildMutationPlot(WaterfallData.out, plotA="not_possible", plotATally="simple", plotALayers=NULL, verbose=FALSE))
})

test_that("buildMutationPlot warns if input to plotATally is unexpected", {

    expect_warning(buildMutationPlot(WaterfallData.out, plotA="frequency", plotATally="not_possible", plotALayers=NULL, verbose=FALSE))
})

test_that("buildMutationPlot warns if plotALayers is not passed as a list", {

    test_plotALayers <- ggplot2::scale_y_log10()
    expect_error(buildMutationPlot(WaterfallData.out, plotA="frequency", plotATally="simple", plotALayers=test_plotALayers, verbose=FALSE))
})

test_that("buildMutationPlot warns if plotALayers does not contain valid ggplot2 layers", {

    test_plotALayers <- list(c("THIS IS A TEST"))
    expect_warning(buildMutationPlot(WaterfallData.out, plotA="frequency", plotATally="simple", plotALayers=test_plotALayers, verbose=FALSE))
})

test_that("buildMutationPlot succesfully adds layers to the plot", {

    test_plotALayers <- list(ggtitle("THIS IS A TEST"), xlab("THIS IS A TEST"), ylab("THIS IS A TEST"))
    buildMutationPlot.out <- buildMutationPlot(WaterfallData.out, plotA="frequency", plotATally="simple", plotALayers=test_plotALayers, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("mutation plot add layers", grid::grid.draw(buildMutationPlot.out))
})

################################# buildGenePlot ################################

context("Waterfall Gene Plot")

test_that("buildGenePlot draws a simple proportion plot correctly", {

    buildGenePlot.out <- buildGenePlot(WaterfallData.out, plotB="proportion", plotBTally="simple", plotBLayers=NULL, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("gene plot proportion simple", grid::grid.draw(buildGenePlot.out))

})

test_that("buildGenePlot draws a complex proportion plot correctly", {

    buildGenePlot.out <- buildGenePlot(WaterfallData.out, plotB="proportion", plotBTally="complex", plotBLayers=NULL, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("gene plot proportion complex", grid::grid.draw(buildGenePlot.out))

})

test_that("buildGenePlot draws a simple frequency plot correctly", {


    buildGenePlot.out <- buildGenePlot(WaterfallData.out, plotB="frequency", plotBTally="simple", plotBLayers=NULL, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("gene plot frequency simple", grid::grid.draw(buildGenePlot.out))

})

test_that("buildGenePlot draws a complex frequency plot correctly", {

    buildGenePlot.out <- buildGenePlot(WaterfallData.out, plotB="frequency", plotBTally="complex", plotBLayers=NULL, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("gene plot frequency complex", grid::grid.draw(buildGenePlot.out))

})

test_that("buildGenePlot works in verbose mode", {

    expect_message(buildGenePlot(WaterfallData.out, plotB="frequency", plotBTally="complex", plotBLayers=NULL, verbose=TRUE))
})

test_that("buildGenePlot warns if input to plotB is unexpected", {

    expect_warning(buildGenePlot(WaterfallData.out, plotB="not_possible", plotBTally="simple", plotBLayers=NULL, verbose=FALSE))
})

test_that("buildGenePlot warns if input to plotBTally is unexpected", {

    expect_warning(buildGenePlot(WaterfallData.out, plotB="frequency", plotBTally="not_possible", plotBLayers=NULL, verbose=FALSE))
})

test_that("buildGenePlot warns if plotBLayers is not passed as a list", {

    test_plotBLayers <- ggplot2::scale_y_log10()
    expect_error(buildGenePlot(WaterfallData.out, plotB="frequency", plotBTally="simple", plotBLayers=test_plotBLayers, verbose=FALSE))
})

test_that("buildGenePlot warns if plotBLayers does not contain valid ggplot2 layers", {

    test_plotBLayers <- list(c("THIS IS A TEST"))
    expect_warning(buildGenePlot(WaterfallData.out, plotB="frequency", plotBTally="simple", plotBLayers=test_plotBLayers, verbose=FALSE))
})

test_that("buildGenePlot succesfully adds layers to the plot", {

    test_plotBLayers <- list(ggtitle("THIS IS A TEST"), xlab("THIS IS A TEST"), ylab("THIS IS A TEST"))
    buildGenePlot.out <- buildGenePlot(WaterfallData.out, plotB="frequency", plotBTally="simple", plotBLayers=test_plotBLayers, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("gene plot add layers", grid::grid.draw(buildGenePlot.out))
})

########################### buildWaterfallPlot #################################

context("Waterfall main Plot")

test_that("buildWaterfallPlot draws a base plot", {

    buildWaterfallPlot.out <- buildWaterfallPlot(WaterfallData.out, gridOverlay =FALSE,
                                                 drop=FALSE, labelSize=5, labelAngle=0,
                                                 sampleNames=TRUE, xTitle=TRUE, plotCLayers=NULL,
                                                 verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("main waterfall plot base", grid::grid.draw(buildWaterfallPlot.out))

})

test_that("buildWaterfallPlot draws a grid", {

    buildWaterfallPlot.out <- buildWaterfallPlot(WaterfallData.out, gridOverlay =TRUE,
                                                 drop=FALSE, labelSize=5, labelAngle=0,
                                                 sampleNames=TRUE, xTitle=TRUE, plotCLayers=NULL,
                                                 verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("main waterfall plot grid true", grid::grid.draw(buildWaterfallPlot.out))

})

test_that("buildWaterfallPlot drops mutations", {

    buildWaterfallPlot.out <- buildWaterfallPlot(WaterfallData.out, gridOverlay =FALSE,
                                                 drop=TRUE, labelSize=5, labelAngle=0,
                                                 sampleNames=TRUE, xTitle=TRUE, plotCLayers=NULL,
                                                 verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("main waterfall plot drop true", grid::grid.draw(buildWaterfallPlot.out))

})

test_that("buildWaterfallPlot doesn't plots samples", {

    buildWaterfallPlot.out <- buildWaterfallPlot(WaterfallData.out, gridOverlay =FALSE,
                                                 drop=FALSE, labelSize=5, labelAngle=0,
                                                 sampleNames=FALSE, xTitle=TRUE, plotCLayers=NULL,
                                                 verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("main waterfall plot sampleNames false", grid::grid.draw(buildWaterfallPlot.out))

})

test_that("buildWaterfallPlot doesn't plot an x-axis title", {

    buildWaterfallPlot.out <- buildWaterfallPlot(WaterfallData.out, gridOverlay =FALSE,
                                                 drop=FALSE, labelSize=5, labelAngle=0,
                                                 sampleNames=FALSE, xTitle=FALSE, plotCLayers=NULL,
                                                 verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("main waterfall plot xtitle false", grid::grid.draw(buildWaterfallPlot.out))

})

test_that("buildWaterfallPlot labels plot cells", {

    # test labeling feature
    WaterfallData.out <- WaterfallData(mafObject, labelColumn="Hugo_Symbol", samples=NULL,
                                       mutationHierarchy=NULL, coverage=44000000, mutation=NULL,
                                       genes=NULL, recurrence=NULL, geneOrder=NULL, geneMax=15,
                                       sampleOrder=NULL, verbose=FALSE)
    buildWaterfallPlot.out <- buildWaterfallPlot(WaterfallData.out, gridOverlay =FALSE,
                                                 drop=FALSE, labelSize=5, labelAngle=0,
                                                 sampleNames=FALSE, xTitle=FALSE, plotCLayers=NULL,
                                                 verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("main waterfall plot label", grid::grid.draw(buildWaterfallPlot.out))

})

test_that("buildWaterfallPlot works in verbose mode", {

    expect_message(buildWaterfallPlot(WaterfallData.out, gridOverlay =FALSE,
                                      drop=FALSE, labelSize=5, labelAngle=0,
                                      sampleNames=TRUE, xTitle=TRUE, plotCLayers=NULL,
                                      verbose=TRUE))
})

test_that("buildWaterfallPlot warns if plotCLayers is not a list", {

    test_plotCLayers <- ggplot2::scale_x_discrete()
    expect_error(buildWaterfallPlot(WaterfallData.out, gridOverlay =FALSE,
                                      drop=FALSE, labelSize=5, labelAngle=0,
                                      sampleNames=TRUE, xTitle=TRUE, plotCLayers=test_plotCLayers,
                                      verbose=TRUE))
})

test_that("buildWaterfallPlot warns if plotCLayers contains an invalid ggplot object", {

    test_plotCLayers <- list(c("THIS IS A TEST"))
    expect_warning(buildWaterfallPlot(WaterfallData.out, gridOverlay =FALSE,
                                      drop=FALSE, labelSize=5, labelAngle=0,
                                      sampleNames=TRUE, xTitle=TRUE, plotCLayers=test_plotCLayers,
                                      verbose=TRUE))
})

test_that("buildWaterfallPlot successfully adds layers to a plot", {

    test_plotCLayers <- list(ggtitle("THIS IS A TEST"), xlab("THIS IS A TEST"), ylab("THIS IS A TEST"))
    buildWaterfallPlot.out <- buildWaterfallPlot(WaterfallData.out, gridOverlay =FALSE,
                                                 drop=FALSE, labelSize=5, labelAngle=0,
                                                 sampleNames=TRUE, xTitle=TRUE, plotCLayers=test_plotCLayers,
                                                 verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("main waterfall add layers", grid::grid.draw(buildWaterfallPlot.out))
})


########################### formatClinicalData #################################

context("Waterfall Clinical Plot setup")

# create simple ClinicalObject for testing
library(ggplot2)
clinData <- data.table::data.table("sample"=c(as.character(unique(getSample(mafObject)$Tumor_Sample_Barcode))), "variable"="a", "value"="b")
clinObject <- Clinical(inputData=clinData, inputFormat = "long", clinicalLayers = theme(axis.text.x=element_text(angle=20)), verbose=FALSE)

formatClinicalData.out <- formatClinicalData(WaterfallData.out, clinical=clinObject, verbose=FALSE)

test_that("formatClinicalData adjusts the Clinical object to match samples in WaterfallData", {
    expected <- levels(getData(WaterfallData.out, name="primaryData")$sample)
    actual <- levels(formatClinicalData.out$sample)
    expect_true(all(expected == actual))
})

test_that("formatClinicalData removes samples not in the WaterfallData", {

    # create simple test
    clinData <- data.table::data.table("sample"=c(as.character(unique(getSample(mafObject)$Tumor_Sample_Barcode)), "test"), "variable"="a", "value"="b")
    clinObject <- Clinical(inputData=clinData, inputFormat = "long", clinicalLayers = theme(axis.text.x=element_text(angle=20)), verbose=FALSE)

    # expect warning
    expect_warning(formatClinicalData(WaterfallData.out, clinical=clinObject, verbose=FALSE), "Removed")

    # expect levels match
    formatClinicalData.out <- suppressWarnings(formatClinicalData(WaterfallData.out, clinical=clinObject, verbose=FALSE))
    expected <- levels(getData(WaterfallData.out, name="primaryData")$sample)
    actual <- levels(formatClinicalData.out$sample)
    expect_true(all(expected == actual))

})

test_that("formatClinicalData fills missing samples not in the Clinical object", {

    # create simple test
    clinData <- data.table::data.table("sample"=c(as.character(unique(getSample(mafObject)$Tumor_Sample_Barcode)[-1])), "variable"="a", "value"="b")
    clinObject <- Clinical(inputData=clinData, inputFormat = "long", clinicalLayers = theme(axis.text.x=element_text(angle=20)), verbose=FALSE)

    # expect warning
    expect_warning(formatClinicalData(WaterfallData.out, clinical=clinObject, verbose=FALSE), "Added")

    # expect levels match
    formatClinicalData.out <- suppressWarnings(formatClinicalData(WaterfallData.out, clinical=clinObject, verbose=FALSE))
    expected <- levels(getData(WaterfallData.out, name="primaryData")$sample)
    actual <- levels(formatClinicalData.out$sample)
    expect_true(all(expected == actual))

})

test_that("formatClinicalData works in verbose mode", {
    expect_message(formatClinicalData(WaterfallData.out, clinical=clinObject, verbose=TRUE))
})

#################### test WaterfallPlot class construction #####################

context("WaterfallPlots Constructor")

WaterfallPlots.out <- WaterfallPlots(WaterfallData.out, clinical=NULL, plotA=NULL,
                                     plotATally="simple", plotALayers=NULL,
                                     plotB=NULL, plotBTally="simple", plotBLayers=NULL,
                                     gridOverlay=FALSE, drop=TRUE, labelSize=5,
                                     labelAngle=0, sampleNames=FALSE, plotCLayers=NULL,
                                     verbose=FALSE)

test_that("WaterfallPlots constructor outputs a S4 class object", {
    expect_s4_class(WaterfallPlots.out, "WaterfallPlots")
})

################################################################################
################ test Waterfall constructor and associated functions ###########

######################### arrangeWaterfallPlot #################################

context("Waterfall Final Plot")

test_that("arrangeWaterfallPlot draws base plot", {

    arrangeWaterfallPlot.out <- arrangeWaterfallPlot(WaterfallPlots.out, sectionHeights=NULL, sectionWidths=NULL, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("final waterfall base", grid::grid.draw(arrangeWaterfallPlot.out))

})

test_that("arrangeWaterfallPlot arranges a top sub-plot", {

    WaterfallPlots.out <- WaterfallPlots(WaterfallData.out, clinical=NULL, plotA="frequency",
                                         plotATally="simple", plotALayers=NULL,
                                         plotB=NULL, plotBTally="simple", plotBLayers=NULL,
                                         gridOverlay=FALSE, drop=TRUE, labelSize=5,
                                         labelAngle=0, sampleNames=FALSE, plotCLayers=NULL,
                                         verbose=FALSE)
    arrangeWaterfallPlot.out <- arrangeWaterfallPlot(WaterfallPlots.out, sectionHeights=NULL, sectionWidths=NULL, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("final waterfall and top sub plot", grid::grid.draw(arrangeWaterfallPlot.out))

})

test_that("arrangeWaterfallPlot arranges a left sub-plot", {

    WaterfallPlots.out <- WaterfallPlots(WaterfallData.out, clinical=NULL, plotA=NULL,
                                         plotATally="simple", plotALayers=NULL,
                                         plotB="frequency", plotBTally="simple", plotBLayers=NULL,
                                         gridOverlay=FALSE, drop=TRUE, labelSize=5,
                                         labelAngle=0, sampleNames=FALSE, plotCLayers=NULL,
                                         verbose=FALSE)
    arrangeWaterfallPlot.out <- arrangeWaterfallPlot(WaterfallPlots.out, sectionHeights=NULL, sectionWidths=NULL, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("final waterfall and left sub plot", grid::grid.draw(arrangeWaterfallPlot.out))

})

test_that("arrangeWaterfallPlot draw both a top and left sub-plot", {

    WaterfallPlots.out <- WaterfallPlots(WaterfallData.out, clinical=NULL, plotA="frequency",
                                         plotATally="simple", plotALayers=NULL,
                                         plotB="frequency", plotBTally="simple", plotBLayers=NULL,
                                         gridOverlay=FALSE, drop=TRUE, labelSize=5,
                                         labelAngle=0, sampleNames=FALSE, plotCLayers=NULL,
                                         verbose=FALSE)
    arrangeWaterfallPlot.out <- arrangeWaterfallPlot(WaterfallPlots.out, sectionHeights=NULL, sectionWidths=NULL, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("final waterfall and top and left sub plot", grid::grid.draw(arrangeWaterfallPlot.out))

})

test_that("arrangeWaterfallPlot draws a clinical sub-plot", {

    WaterfallPlots.out <- WaterfallPlots(WaterfallData.out, clinical=clinObject, plotA="frequency",
                                         plotATally="simple", plotALayers=NULL,
                                         plotB="frequency", plotBTally="simple", plotBLayers=NULL,
                                         gridOverlay=FALSE, drop=TRUE, labelSize=5,
                                         labelAngle=0, sampleNames=FALSE, plotCLayers=NULL,
                                         verbose=FALSE)

    arrangeWaterfallPlot.out <- arrangeWaterfallPlot(WaterfallPlots.out, sectionHeights=NULL, sectionWidths=NULL, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("final waterfall and top and left and bottom sub plot", grid::grid.draw(arrangeWaterfallPlot.out))

})

test_that("arrangeWaterfallPlot re-sizes plots based on sectionWidths", {

    WaterfallPlots.out <- WaterfallPlots(WaterfallData.out, clinical=NULL, plotA=NULL,
                                         plotATally="simple", plotALayers=NULL,
                                         plotB="frequency", plotBTally="simple", plotBLayers=NULL,
                                         gridOverlay=FALSE, drop=TRUE, labelSize=5,
                                         labelAngle=0, sampleNames=FALSE, plotCLayers=NULL,
                                         verbose=FALSE)
    arrangeWaterfallPlot.out <- arrangeWaterfallPlot(WaterfallPlots.out, sectionHeights=NULL, sectionWidths=c(.5, .5), verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("final waterfall and left sub plot alter sectionWidth", grid::grid.draw(arrangeWaterfallPlot.out))

})

test_that("arrangeWaterfallPlot re-sizes plots based on sectionHeights", {

    WaterfallPlots.out <- WaterfallPlots(WaterfallData.out, clinical=NULL, plotA="frequency",
                                         plotATally="simple", plotALayers=NULL,
                                         plotB=NULL, plotBTally="simple", plotBLayers=NULL,
                                         gridOverlay=FALSE, drop=TRUE, labelSize=5,
                                         labelAngle=0, sampleNames=FALSE, plotCLayers=NULL,
                                         verbose=FALSE)
    arrangeWaterfallPlot.out <- arrangeWaterfallPlot(WaterfallPlots.out, sectionHeights=c(.5, .5), sectionWidths=NULL, verbose=FALSE)
    skip_on_bioc()
    vdiffr::expect_doppelganger("final waterfall and top sub plot alter sectionHeights", grid::grid.draw(arrangeWaterfallPlot.out))


    expect_warning(arrangeWaterfallPlot(WaterfallPlots.out, sectionHeights=c(.5, .5, .5), sectionWidths=NULL, verbose=FALSE))
})

################################################################################
############## Test Waterfall Constructor and accessors ########################

context("Waterfall Constructor")

Waterfall.out <- Waterfall(mafObject, labelColumn=NULL, samples=NULL, coverage=NULL,
                           mutation=NULL, genes=NULL, mutationHierarchy=NULL,
                           recurrence=NULL, geneOrder=NULL, geneMax=5,
                           sampleOrder=NULL, plotA=c("frequency", "burden", NULL),
                           plotATally=c("simple", "complex"), plotALayers=NULL,
                           plotB=c("proportion", "frequency", NULL),
                           plotBTally=c("simple", "complex"), plotBLayers=NULL,
                           gridOverlay=FALSE, drop=TRUE, labelSize=5, labelAngle=0,
                           sampleNames=TRUE, clinical=NULL, sectionHeights=NULL,
                           sectionWidths=NULL, verbose=FALSE, plotCLayers=NULL)

test_that("Waterfall constructor outputs a S4 class object", {
    expect_s4_class(Waterfall.out, "Waterfall")
})

test_that("drawPlot constructs a waterfall plot from grob objects in Waterfall object", {

    skip_on_bioc()
    vdiffr::expect_doppelganger("drawPlot waterfall", drawPlot(Waterfall.out))
})

test_that("Waterfall adds genes into the plot which are not in the data", {

        Waterfall.out <- suppressWarnings(Waterfall(mafObject, labelColumn=NULL, samples=NULL, coverage=NULL,
                               mutation=NULL, genes=c("not_here_1", "not_here_2"),
                               mutationHierarchy=NULL,
                               recurrence=.4, geneOrder=NULL, geneMax=NULL,
                               sampleOrder=NULL, plotA=c("frequency", "burden", NULL),
                               plotATally=c("simple", "complex"), plotALayers=NULL,
                               plotB=c("proportion", "frequency", NULL),
                               plotBTally=c("simple", "complex"), plotBLayers=NULL,
                               gridOverlay=FALSE, drop=TRUE, labelSize=5, labelAngle=0,
                               sampleNames=TRUE, clinical=NULL, sectionHeights=NULL,
                               sectionWidths=NULL, verbose=FALSE, plotCLayers=NULL))
    skip_on_bioc()
    vdiffr::expect_doppelganger("addgene waterfall", drawPlot(Waterfall.out))
})

########################### getData ############################################

context("Waterfall Accessors")

test_that("getData outputs error if no name or index is given", {

    expect_error(getData(Waterfall.out))

})

test_that("getData outputs error if index exceeds the number of slots", {

    expect_error(getData(Waterfall.out, index=10))

})

test_that("getData outputs error if supplied name is not a valid slot name", {

    expect_error(getData(Waterfall.out, name="shouldNotexist"))

})

test_that("getData retrieves specified slot data correctly", {

    expect_s3_class(getData(Waterfall.out, index=1), "data.table")
    expect_equivalent(getData(Waterfall.out, name="primaryData"), getData(Waterfall.out, index=1))

    expect_s3_class(getData(Waterfall.out, index=2), "data.table")
    expect_equivalent(getData(Waterfall.out, name="simpleMutationCounts"), getData(Waterfall.out, index=2))

    expect_s3_class(getData(Waterfall.out, index=3), "data.table")
    expect_equivalent(getData(Waterfall.out, name="complexMutationCounts"), getData(Waterfall.out, index=3))

    expect_s3_class(getData(Waterfall.out, index=4), "data.table")
    expect_equivalent(getData(Waterfall.out, name="geneData"), getData(Waterfall.out, index=4))

    expect_s3_class(getData(Waterfall.out, index=5), "data.table")
    expect_equivalent(getData(Waterfall.out, name="mutationHierarchy"), getData(Waterfall.out, index=5))

})

################## getGrob #####################################################

test_that("getGrob outputs error if index is out of bounds", {

    expect_error(getGrob(Waterfall.out, index=10))
})

test_that("getGrob successfully retrieves grob objects from Waterfall object", {

    expect_s3_class(getGrob(Waterfall.out, index=1), "gtable")
    expect_s3_class(getGrob(Waterfall.out, index=2), "gtable")
    expect_s3_class(getGrob(Waterfall.out, index=3), "gtable")
    expect_s3_class(getGrob(Waterfall.out, index=4), "gtable")
})
