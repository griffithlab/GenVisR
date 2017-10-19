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



############################## recurrenceSubset ################################
