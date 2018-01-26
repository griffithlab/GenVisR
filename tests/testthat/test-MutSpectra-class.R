# packages needed
library(ggplot2)

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

test_that("annoMutSpectra works in verbose mode", {
    expect_message(annoMutSpectra(toMutSpectra.out, verbose=TRUE))
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

test_that("calcMutSpectra works in verbose mode", {
    expect_message(calcMutSpectra(annoMutSpectra.out, verbose=TRUE))
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

test_that("sortSamples warns if incorrect input is detected to sorting", {
    expect_warning(sortSamples(calcMutSpectra.out, sorting=as.factor("sample"), verbose=FALSE))
})

test_that("sortSamples works in verbose mode", {
    sampleOrder <- c("FLX005-Naive", "FLX007-Naive", "FLX001-Naive", "FLX004-Naive", "FLX003-Naive")
    expect_message(sortSamples(calcMutSpectra.out, sorting="sample", verbose=TRUE))
    expect_message(sortSamples(calcMutSpectra.out, sorting="mutation", verbose=TRUE))
    expect_message(sortSamples(calcMutSpectra.out, sorting=sampleOrder, verbose=TRUE))
})

#### test MutSpectraPrimaryData class construction with various parameters #####

MutSpectraPrimaryData.out <- MutSpectraPrimaryData(gmsObject, BSgenome=NULL, sorting=NULL, verbose=FALSE)

test_that("MutSpectraPrimaryData outputs a S4 class object", {
    expect_s4_class(MutSpectraPrimaryData.out, "MutSpectraPrimaryData")
})

################################################################################
######### test the MutSpectraPlots constructor and various methods #############
################################################################################

################## test formatClinicalData #####################################

context("MutSpectra Clinical Plot setup")

# create simple ClinicalObject for testing
library(ggplot2)
clinData <- data.table::data.table("sample"=c(as.character(unique(getSample(gmsObject)$sample))), "variable"="a", "value"="b")
clinObject <- Clinical(inputData=clinData, inputFormat = "long", clinicalLayers = theme(axis.text.x=element_text(angle=20)), verbose=FALSE)

formatClinicalData.out <- formatClinicalData(MutSpectraPrimaryData.out, clinical=clinObject, verbose=FALSE)

test_that("formatClinicalData adjusts the Clinical object to match samples in MutSpectraPrimaryData", {
    expected <- levels(getData(MutSpectraPrimaryData.out, name="primaryData")$sample)
    actual <- levels(formatClinicalData.out$sample)
    expect_true(all(expected == actual))
})

test_that("formatClinicalData removes samples not in the MutSpectraPrimaryData", {
    
    # create simple test
    clinData <- data.table::data.table("sample"=c(as.character(unique(getSample(gmsObject)$sample)), "test"), "variable"="a", "value"="b")
    clinObject <- Clinical(inputData=clinData, inputFormat = "long", clinicalLayers = theme(axis.text.x=element_text(angle=20)), verbose=FALSE)
    
    # expect warning
    expect_warning(formatClinicalData(MutSpectraPrimaryData.out, clinical=clinObject, verbose=FALSE), "Removed")
    
    # expect levels match
    formatClinicalData.out <- suppressWarnings(formatClinicalData(MutSpectraPrimaryData.out, clinical=clinObject, verbose=FALSE))
    expected <- levels(getData(MutSpectraPrimaryData.out, name="primaryData")$sample)
    actual <- levels(formatClinicalData.out$sample)
    expect_true(all(expected == actual))
    
})

test_that("formatClinicalData fills missing samples not in the Clinical object", {
    
    # create simple test
    clinData <- data.table::data.table("sample"=c(as.character(unique(getSample(gmsObject)$sample)[-1])), "variable"="a", "value"="b")
    clinObject <- Clinical(inputData=clinData, inputFormat = "long", clinicalLayers = theme(axis.text.x=element_text(angle=20)), verbose=FALSE)
    
    # expect warning
    expect_warning(formatClinicalData(MutSpectraPrimaryData.out, clinical=clinObject, verbose=FALSE), "Added")
    
    # expect levels match
    formatClinicalData.out <- suppressWarnings(formatClinicalData(MutSpectraPrimaryData.out, clinical=clinObject, verbose=FALSE))
    expected <- levels(getData(MutSpectraPrimaryData.out, name="primaryData")$sample)
    actual <- levels(formatClinicalData.out$sample)
    expect_true(all(expected == actual))
    
})

test_that("formatClinicalData works in verbose mode", {
    expect_message(formatClinicalData(MutSpectraPrimaryData.out, clinical=clinObject, verbose=TRUE))
})

######################### test buildFrequencyPlot ##############################

context("MutSpectra Frequency Plot")


test_that("buildFrequencyPlot constructs a plot based on frequencies", {
    
    skip_on_bioc()
    
    buildFrequencyPlot.out <- buildFrequencyPlot(MutSpectraPrimaryData.out, plotALayers=NULL, palette=NULL, verbose=FALSE)
    vdiffr::expect_doppelganger("mutspectra frequency plot", grid::grid.draw(buildFrequencyPlot.out))
    
})

test_that("buildFrequencyPlot is able to add layers to the plot", {
    
    skip_on_bioc()
    
    test_plotALayers <- list(ggplot2::geom_hline(yintercept=c(30), colour="black", size=2), ggplot2::geom_vline(xintercept=c(2), colour="black", size=2))
    buildFrequencyPlot.out <- buildFrequencyPlot(MutSpectraPrimaryData.out, plotALayers=test_plotALayers, palette=NULL, verbose=FALSE)
    vdiffr::expect_doppelganger("mutspectra frequency plot add layer", grid::grid.draw(buildFrequencyPlot.out))
      
})

test_that("buildFrequencyPlot is able to add a custom pallete to the plot", {
    
    skip_on_bioc()
    
    testPallete <- c("red", "blue", "green", "black", "darkorchid4", "seagreen3")
    buildFrequencyPlot.out <- buildFrequencyPlot(MutSpectraPrimaryData.out, plotALayers=NULL, palette=testPallete, verbose=FALSE)
    vdiffr::expect_doppelganger("mutspectra frequency plot custom pallette", grid::grid.draw(buildFrequencyPlot.out))
})

test_that("buildFrequencyPlot warns if plotALayers is not passed as a list", {
    
    test_plotALayers <- ggplot2::geom_hline(yintercept=c(30), colour="black", size=2)
    expect_error(buildFrequencyPlot(MutSpectraPrimaryData.out, plotALayers=test_plotALayers, palette=NULL, verbose=FALSE))
    
})

test_that("buildFrequencyPlot warns if plotALayers does not contain valid ggplot2 layers", {
    
    test_plotALayers <- list(c("THIS IS A TEST"))
    expect_warning(buildFrequencyPlot(MutSpectraPrimaryData.out, plotALayers=test_plotALayers, palette=NULL, verbose=FALSE))
})

test_that("buildFrequencyPlot works in verbose mode", {
    expect_message(buildFrequencyPlot(MutSpectraPrimaryData.out, plotALayers=NULL, palette=NULL, verbose=TRUE))
})

test_that("buildFrequencyPlot warns if a pallete is not the correct length", {
    testPallete <- c("red", "blue", "green", "black", "darkorchid4", "seagreen3", "salmon")
    expect_warning(buildFrequencyPlot(MutSpectraPrimaryData.out, plotALayers=NULL, palette=testPallete, verbose=FALSE))
})

######################### test buildProportionPlot #############################

context("MutSpectra Proportion Plot")

test_that("buildProportionPlot constructs a plot based on Proportions", {
    
    skip_on_bioc()
    
    buildProportionPlot.out <- buildProportionPlot(MutSpectraPrimaryData.out, sampleNames=TRUE, plotBLayers=NULL, palette=NULL, verbose=FALSE)
    vdiffr::expect_doppelganger("mutspectra proportion plot", grid::grid.draw(buildProportionPlot.out))
    
})

test_that("buildProportionPlot is able to add layers to the plot", {
    
    skip_on_bioc()
    
    test_plotBLayers <- list(ggplot2::geom_hline(yintercept=c(.5), colour="black", size=2), ggplot2::geom_vline(xintercept=c(2), colour="black", size=2))
    buildProportionPlot.out <- buildProportionPlot(MutSpectraPrimaryData.out, sampleNames=TRUE, plotBLayers=test_plotBLayers, palette=NULL, verbose=FALSE)
    vdiffr::expect_doppelganger("mutspectra proportion plot add layer", grid::grid.draw(buildProportionPlot.out))
    
})

test_that("buildProportionPlot is able to add a custom pallete to the plot", {
    
    skip_on_bioc()
    
    testPallete <- c("red", "blue", "green", "black", "darkorchid4", "seagreen3")
    buildProportionPlot.out <- buildProportionPlot(MutSpectraPrimaryData.out, sampleNames=FALSE, plotBLayers=NULL, palette=testPallete, verbose=FALSE)
    vdiffr::expect_doppelganger("mutspectra proportion plot custom pallette", grid::grid.draw(buildProportionPlot.out))
})

test_that("buildProportionPlot warns if plotBLayers is not passed as a list", {
    
    test_plotBLayers <- ggplot2::geom_hline(yintercept=c(.5), colour="black", size=2)
    expect_error(buildProportionPlot(MutSpectraPrimaryData.out, sampleNames=TRUE, plotBLayers=test_plotBLayers, palette=NULL, verbose=FALSE))
    
})

test_that("buildProportionPlot warns if plotBLayers does not contain valid ggplot2 layers", {
    
    test_plotBLayers <- list(c("THIS IS A TEST"))
    expect_warning(buildProportionPlot(MutSpectraPrimaryData.out, sampleNames=TRUE, plotBLayers=test_plotBLayers, palette=NULL, verbose=FALSE))
})

test_that("buildProportionPlot works in verbose mode", {
    expect_message(buildProportionPlot(MutSpectraPrimaryData.out, sampleNames=FALSE, plotBLayers=NULL, palette=NULL, verbose=TRUE))
})

test_that("buildProportionPlot warns if a pallete is not the correct length", {
    testPallete <- c("red", "blue", "green", "black", "darkorchid4", "seagreen3", "salmon")
    expect_warning(buildProportionPlot(MutSpectraPrimaryData.out, sampleNames=FALSE, plotBLayers=NULL, palette=testPallete, verbose=FALSE))
})


######################### test MutSpectraPlots class construction ##############

context("MutSpectraPlots Constructor")

MutSpectraPlots.out <- MutSpectraPlots(MutSpectraPrimaryData.out, clinical=NULL, plotALayers=NULL,
                                       plotBLayers=NULL, plotCLayers=NULL, sampleNames=FALSE,
                                       palette=NULL, verbose=FALSE)

test_that("MutSpectraPlots constructor outputs a s4 class object", {
    
    expect_s4_class(MutSpectraPlots.out, "MutSpectraPlots")
    
})

################################################################################
############# test MutSpectra constructor and associated functions #############

context("MutSpectra Final Plot")

################# arrangeMutSpectraPlot ########################################

test_that("arrangeMutSpectraPlot plots a base plot", {
    
    skip_on_bioc()
    
    arrangeMutSpectraPlot.out <- arrangeMutSpectraPlot(MutSpectraPlots.out, sectionHeights=NULL, verbose=FALSE)
    vdiffr::expect_doppelganger("final MutSpectra base", grid::grid.draw(arrangeMutSpectraPlot.out))
    
})

test_that("arrangeMutSpectraPlots alters section heights", {
    
    skip_on_bioc()
    
    arrangeMutSpectraPlot.out <- arrangeMutSpectraPlot(MutSpectraPlots.out, sectionHeights=c(1, 3), verbose=FALSE)
    vdiffr::expect_doppelganger("final MutSpectra alter section heights", grid::grid.draw(arrangeMutSpectraPlot.out))
    
})

test_that("arrangeMutSpectraPlots correctly warns if section heights does not match the number of plots", {
    expect_warning(arrangeMutSpectraPlot(MutSpectraPlots.out, sectionHeights=c(1), verbose=FALSE))
})

test_that("arrangeMutSpectraPlots correctly warns if section heights is not numeric", {
    expect_warning(arrangeMutSpectraPlot(MutSpectraPlots.out, sectionHeights=c("A"), verbose=FALSE))
})

################################################################################
############## Test MutSpectra Constructor and accessors ########################

MutSpectra.out <- MutSpectra(gmsObject, BSgenome=NULL, sorting=NULL, palette=NULL,
                             clinical=NULL, sectionHeights=NULL, sampleNames=TRUE,
                             verbose=FALSE, plotALayers=NULL, plotBLayers=NULL,
                             plotCLayers=NULL)

test_that("MutSpectra constructor outputs a S4 class object", {
    expect_s4_class(MutSpectra.out, "MutSpectra")
})

test_that("drawPlot constructs a MutSpectra plot from grob objects in MutSpectra object", {
    
    skip_on_bioc()
    
    vdiffr::expect_doppelganger("drawPlot MutSpectra", drawPlot(MutSpectra.out))
})

test_that("MutSpectra correctly adds clinical data", {
    
    skip_on_bioc()

    MutSpectra.out <- MutSpectra(gmsObject, BSgenome=NULL, sorting=NULL, palette=NULL,
                                 clinical=clinObject, sectionHeights=NULL, sampleNames=TRUE,
                                 verbose=FALSE, plotALayers=NULL, plotBLayers=NULL,
                                 plotCLayers=NULL)
    vdiffr::expect_doppelganger("MutSpectra Clinical", drawPlot(MutSpectra.out))

})

################## getGrob #####################################################

context("MutSpectra accessors")

test_that("getGrob outputs error if index is out of bounds", {

    expect_error(getGrob(MutSpectra.out, index=10))
})

test_that("getGrob successfully retrieves grob objects from MutSpectra object", {

    expect_s3_class(getGrob(MutSpectra.out, index=1), "gtable")
    expect_s3_class(getGrob(MutSpectra.out, index=2), "gtable")
    expect_s3_class(getGrob(MutSpectra.out, index=3), "gtable")
    expect_s3_class(getGrob(MutSpectra.out, index=4), "gtable")
})

################## getData #####################################################

test_that("getData outputs error if no name or index is given", {
    
    expect_error(getData(MutSpectra.out))
    
})

test_that("getData outputs error if index exceeds the number of slots", {
    
    expect_error(getData(MutSpectra.out, index=10))
    
})

test_that("getData outputs error if supplied name is not a valid slot name", {
    
    expect_error(getData(MutSpectra.out, name="shouldNotexist"))
    
})

test_that("getData retrieves specified slot data correctly", {
    
    expect_s3_class(getData(MutSpectra.out, index=1), "data.table")
    expect_equivalent(getData(MutSpectra.out, name="primaryData"), getData(MutSpectra.out, index=1))
    
    expect_s3_class(getData(MutSpectra.out, index=2), "data.table")
    expect_equivalent(getData(MutSpectra.out, name="ClinicalData"), getData(MutSpectra.out, index=2))
    
})