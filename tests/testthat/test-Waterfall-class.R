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

############################### rmvMutation ####################################

test_that("rmvMutation warns if a mutation was specified to be kept but is not found", {
    expect_warning(rmvMutation(toWaterfall.out, mutation=c("invalid"), verbose=FALSE), "were not found")
})

rmvMutation.out <- rmvMutation(toWaterfall.out, mutation=c("Missense_Mutation"), verbose=FALSE)
test_that("rmvMutation keeps only mutation specified to be kept", {
    expect_true(all(rmvMutation.out$mutation == "Missense_Mutation"))
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

############################## geneSubset ######################################

geneSubset.out <- geneSubset(mutHierarchySubset.out, genes=NULL, verbose=FALSE)
test_that("geneSubset keeps samples by outputing an NA value in addition to genes", {
    
    expect_true(any(is.na(geneSubset.out)))
    
    geneSubset.out <- geneSubset(mutHierarchySubset.out, genes=c("ARC"), verbose=FALSE)
    expect_true(any(is.na(geneSubset.out)))
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

######################### geneFilter ###########################################

keepGenes <- unique(c(geneSubset.out, recurrenceSubset.out))
geneFilter.out <- geneFilter(mutHierarchySubset.out, genes=keepGenes, verbose=FALSE)

test_that("geneFilter correctly subsets data to only genes specified", {
    expect_true(all(geneFilter.out$gene %in% keepGenes))
})

test_that("geneFilter correctly identifies if a gene is specified to be filtered but is not found", {
    keepGenes <- c(keepGenes, "not_in_data")
    expect_warning(geneFilter(mutHierarchySubset.out, genes=keepGenes, verbose=FALSE), "not found in")
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

######### test WaterfallData class construction with various parameters ########

# TODO figure out a better way to do or drop it
messages <- capture.output(WaterfallData(mafObject, labelColumn = NULL, samples = NULL,
                                   mutationHierarchy = NULL, coverage = NULL, mutation = NULL,
                                   genes = NULL, recurrence = NULL, geneOrder = NULL, geneMax = NULL,
                                   sampleOrder = NULL, verbose = TRUE), type="message")

test_that("WaterfallData constructor outputs a minimum number of status messages", {
    expect_true(length(messages) > 5)
})
