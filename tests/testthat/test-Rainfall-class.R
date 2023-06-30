# packges needed
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome <- BSgenome.Hsapiens.UCSC.hg19

# get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")
testFile <- Sys.glob(paste0(testFileDir, "/*.vep"))

# define the object for testing
vepObject <- VEP(testFile)

################################################################################
###### Test Rainfall class and associated functions in constructor #############
################################################################################

context("Rainfall Constructor")

###################### toRainfall ##############################################

toRainfall.out <- suppressWarnings(toRainfall(vepObject, BSgenome=BSgenome, verbose=FALSE))

###################### test method for datatable and dataframe #################

test_that("toRainfall removes duplicate genomic mutations", {
    
    datatableObject <- data.table::data.table("sample"=rep("test", 3), "chromosome"=rep(1, 3),
                                              "start"=c(rep(12345, 2), 112358), "stop"=c(rep(12345, 2), 112358),
                                              "refAllele"=rep("A", 3), "variantAllele"=rep("G", 3))
    toRainfall.out <- suppressWarnings(toRainfall(datatableObject, BSgenome=NULL, verbose=FALSE))
    
    expect <- 2
    actual <- nrow(toRainfall.out)
    expect_equal(expect, actual)
})

test_that("toRainfall removes entries with no mutation", {
    
    datatableObject <- data.table::data.table("sample"=rep("test", 3), "chromosome"=rep(1, 3),
                                              "start"=c(rep(12345, 2), 112358), "stop"=c(rep(12345, 2), 112358),
                                              "refAllele"=rep("A", 3), "variantAllele"=c(rep("A", 2), "G"))
    toRainfall.out <- suppressWarnings(toRainfall(datatableObject, BSgenome=NULL, verbose=FALSE))
    
    expect <- 1
    actual <- nrow(toRainfall.out)
    expect_equal(expect, actual)
})

test_that("toRainfall errors if the proper columns are not present", {
    
    datatableObject <- data.table::data.table("wrong"=rep("test", 3), "chromosome"=rep(1, 3),
                                              "start"=c(rep(12345, 2), 112358), "stop"=c(rep(12345, 2), 112358),
                                              "refAllele"=rep("A", 3), "variantAllele"=c(rep("A", 2), "G"))
    expect_error(toRainfall(datatableObject, BSgenome=NULL, verbose=FALSE), "sample")
})

test_that("toRainfall works in verbose mode", {
    
    datatableObject <- data.table::data.table("sample"=rep("test", 2), "chromosome"=rep(1, 2),
                                              "start"=c(12345, 112358), "stop"=c(12345, 112358),
                                              "refAllele"=rep("A", 2), "variantAllele"=rep("G", 2))
    expect_message(toRainfall(datatableObject, BSgenome=NULL, verbose=TRUE))
    
    dataframeObject <- as.data.frame(datatableObject)
    expect_message(toRainfall(dataframeObject, BSgenome=NULL, verbose=TRUE), "converting")
})

test_that("toRainfall returns object of the proper type", {
    
    datatableObject <- data.table::data.table("sample"=rep("test", 2), "chromosome"=rep(1, 2),
                                              "start"=c(12345, 112358), "stop"=c(12345, 112358),
                                              "refAllele"=rep("A", 2), "variantAllele"=rep("G", 2))
    expect_s3_class(toRainfall(datatableObject, BSgenome=NULL, verbose=TRUE), "data.table")
    
    dataframeObject <- as.data.frame(datatableObject)
    expect_s3_class(toRainfall(dataframeObject, BSgenome=NULL, verbose=TRUE), "data.table")
})

####################### annoRainfall ###########################################
annoRainfall.out <- annoRainfall(toRainfall.out, verbose=FALSE)

test_that("annoRainfall correctly annotates transitions and transversions", {

    actual <- as.character(unlist(annoRainfall.out[annoRainfall.out$refAllele == "A" & annoRainfall.out$variantAllele == "G",][1,"trans_tranv"]))
    expected <- "A->G or T->C (TI)"
    expect_equal(actual, expected)

    actual <- as.character(unlist(annoRainfall.out[annoRainfall.out$refAllele == "A" & annoRainfall.out$variantAllele == "C",][1,"trans_tranv"]))
    expected <- "A->C or T->G (TV)"
    expect_equal(actual, expected)

    actual <- as.character(unlist(annoRainfall.out[annoRainfall.out$refAllele == "T" & annoRainfall.out$variantAllele == "G",][1,"trans_tranv"]))
    expected <- "A->C or T->G (TV)"
    expect_equal(actual, expected)

})

test_that("annoRainfall works in verbose mode", {

    expect_message(annoRainfall(toRainfall.out, verbose=TRUE))
})

# ####################### calcMutDist ############################################
calcMutDist.out <- calcMutDist(annoRainfall.out, verbose=FALSE)

test_that("calcMutDist distance calculations reset when jumping between chromosomes", {

    expected <- c(NA, 7.067636, 6.287262)
    actual <- calcMutDist.out[calcMutDist.out$chromosome == "chr20",]
    actual <- actual[order(actual$stop),]
    actual <- actual$log10_difference
    expect_equivalent(expected, actual, tolerance=1e-2)
})

test_that("calcMutDist can calculate distances for duplicate genomic locations", {

    annoRainfall.out <- annoRainfall.out[annoRainfall.out$chromosome == "chr20",]
    annoRainfall.out <- rbind(annoRainfall.out, annoRainfall.out[3,])
    actual <- calcMutDist(annoRainfall.out, verbose=FALSE)
    actual <- actual[actual$start == 60921746]$log10_difference[2]
    expected <- 0
    expect_equivalent(actual, expected)
})

test_that("calcMutDist works in verbose mode", {
    expect_message(calcMutDist(annoRainfall.out, verbose=TRUE))
})

# ####################### chrSubset ##############################################
chrSubset.out <- chrSubset(calcMutDist.out, chromosomes=c("chr8"), verbose=FALSE)

test_that("chrSubset removes all chromosomes not specified", {
    
    actual <- length(unique(chrSubset.out$chromosome))
    expected <- 1
    
    expect_equal(expected, actual)
    
    actual <- chrSubset(calcMutDist.out, chromosomes=c("chr8", "chr1"), verbose=FALSE)
    actual <- length(unique(actual$chromosome))
    expected <- 2
    
    expect_equal(expected, actual)
})

test_that("chrSubset warns if a given chromosome is not found", {
    
    expect_warning(chrSubset(calcMutDist.out, chromosomes=c("chr8", "does_not_exist"), verbose=FALSE))
})

test_that("chrSubset errors if a chromosome subset would eliminate all entries", {
    
    expect_error(expect_warning(chrSubset(calcMutDist.out, chromosomes=c("does_not_exist"), verbose=FALSE)))
})

test_that("chrSubset works in verbose mode", {
    
    expect_message(chrSubset(calcMutDist.out, chromosomes=c("chr8"), verbose=TRUE))
})

####################### annoGenomeCoord ########################################
annoGenomeCoord.out <- annoGenomeCoord(chrSubset.out, BSgenome=BSgenome, verbose=FALSE)

test_that("annoGenomeCoord successfully adds chromosome boundaries from the BSgenome object", {
    
    expected <- 1
    actual <- unique(annoGenomeCoord.out[annoGenomeCoord.out$origin=="chrStart",]$start)
    expect_equal(expected, actual)
    
    actual <- unique(annoGenomeCoord.out[annoGenomeCoord.out$origin=="chrStart",]$stop)
    expect_equal(expected, actual)
    
    # based on chr8
    expected <- 146364022
    actual <- unique(annoGenomeCoord.out[annoGenomeCoord.out$origin=="chrStop",]$start)
    expect_equal(expected, actual)
    
    actual <- unique(annoGenomeCoord.out[annoGenomeCoord.out$origin=="chrStop",]$stop)
    expect_equal(expected, actual)
})

test_that("annoGenomeCoord warns if BSgenome is not set", {
    
    expect_warning(annoGenomeCoord(chrSubset.out, BSgenome=NULL, verbose=FALSE))
})

test_that("annoGenomeCoord checks if appending chr causes a BSgenome match", {
    
    chrSubset.out$chromosome <- gsub("chr", "", chrSubset.out$chromosome)
    expect_warning(annoGenomeCoord(chrSubset.out, BSgenome=BSgenome, verbose=FALSE), "The following")
})

test_that("annoGenomeCoord stops if an unexpected input is given to BSgenome", {
    
    expect_error(annoGenomeCoord(chrSubset.out, BSgenome="character", verbose=FALSE))
})

test_that("annoGenomeCoord works in verbose mode", {
    
    expect_message(annoGenomeCoord(chrSubset.out, BSgenome=BSgenome, verbose=TRUE))
})

######################## formatSample ##########################################
formatSample.out <- formatSample(annoGenomeCoord.out, sample="FLX0070Naive", verbose=FALSE)

test_that("formatSample correctly subsets samples if the sample parameter is set", {
    
    expected <- "FLX0070Naive"
    actual <- as.character(unique(formatSample.out$sample))
    expect_equal(expected, actual)
})

test_that("formatSample checks if sample is a character vector", {
    
    expect_warning(formatSample(annoGenomeCoord.out, sample=as.factor("FLX0070Naive"), verbose=FALSE))
})

test_that("formatSample errors if a sample is given that does not exist", {
    
    expect_error(formatSample(annoGenomeCoord.out, sample="does_not_exist", verbose=FALSE))
})

test_that("formatSample works in verbose mode", {
    
    expect_message(formatSample(annoGenomeCoord.out, sample="FLX0070Naive", verbose=TRUE))
})

################# test RainfallPrimaryData constructor #########################

RainfallPrimaryData.out <- suppressWarnings(RainfallPrimaryData(vepObject, BSgenome=BSgenome, sample=NULL, chromosomes=c("chr1", "chr2") , verbose=FALSE))

test_that("RainfallPrimaryData outputs a S4 class object", {
    
    expect_s4_class(RainfallPrimaryData.out, "RainfallPrimaryData")
})

################################################################################
############## test RainfallPlots constructor and various methods ##############

############# test buildRainfallPlot ###########################################

context("Rainfall Main Plot")

test_that("buildRainfallPlot contructs a plot", {
    
    skip_on_travis()
    
    buildRainfallPlot.out <- buildRainfallPlot(RainfallPrimaryData.out, palette=NULL, pointSize=NULL, plotALayers=NULL, verbose=FALSE)
    skip_on_bioc(vdiffr::expect_doppelganger("Rainfall Main Plot", grid::grid.draw(buildRainfallPlot.out)))
})

test_that("buildRainfallPlot is able to add layers to the plot", {
    
    skip_on_travis()
    
    test_plotALayers <- list(ggplot2::geom_hline(yintercept=c(5), colour="black", size=2), ggplot2::geom_vline(xintercept=c(1.5e8), colour="black", size=2))
    buildRainfallPlot.out <- buildRainfallPlot(RainfallPrimaryData.out, palette=NULL, pointSize=NULL, plotALayers=test_plotALayers , verbose=FALSE)
    skip_on_bioc(vdiffr::expect_doppelganger("Rainfall Plot add layer", grid::grid.draw(buildRainfallPlot.out)))
})

test_that("buildRainfallPlot aesthetic options work", {
    
    skip_on_travis()
   
    colorPalette <- c("red", "blue", "green", "yellow", "orange", "purple", "black")
    buildRainfallPlot.out <- buildRainfallPlot(RainfallPrimaryData.out, palette=colorPalette, pointSize=10, plotALayers=NULL , verbose=FALSE)
    skip_on_bioc(vdiffr::expect_doppelganger("Rainfall Plot aesthetic options", grid::grid.draw(buildRainfallPlot.out)))
})

test_that("buildRainfallPlot warns if plotALayers is not passed as a list", {
    
    test_plotALayers <- ggplot2::geom_hline(yintercept=c(30), colour="black", size=2)
    expect_error(buildRainfallPlot(RainfallPrimaryData.out, palette=NULL, pointSize=NULL, plotALayers=test_plotALayers , verbose=FALSE))
})

test_that("buildRainfallPlot warns if plotALayers does not contain valid ggplot2 layers", {
    
    test_plotALayers <- list(c("THIS IS A TEST"))
    expect_warning(buildRainfallPlot(RainfallPrimaryData.out, palette=NULL, pointSize=NULL, plotALayers=test_plotALayers , verbose=FALSE))
})

test_that("buildRainfallPlot checks aesthetic parameters", {
    colorPalette <- c("red", "blue", "green", "yellow", "orange", "purple", "black", "gray")
    expect_warning(buildRainfallPlot(RainfallPrimaryData.out, palette=NULL, pointSize=as.character(c(1)), plotALayers=NULL, verbose=FALSE))
    expect_warning(buildRainfallPlot(RainfallPrimaryData.out, palette=NULL, pointSize=c(1, 2), plotALayers=NULL, verbose=FALSE))
    expect_warning(buildRainfallPlot(RainfallPrimaryData.out, palette=colorPalette, pointSize=NULL, plotALayers=NULL, verbose=FALSE))
})

test_that("buildRainfallPlot works in verbose mode", {
    expect_message(buildRainfallPlot(RainfallPrimaryData.out, palette=NULL, pointSize=NULL, plotALayers=NULL, verbose=TRUE))
})

############# test buildDensityPlots ###########################################

context("Rainfall Density Plot")

test_that("buildDensityPlot contructs a plot", {
    
    skip_on_travis()
    
    buildDensityPlot.out <- buildDensityPlot(RainfallPrimaryData.out, plotBLayers=NULL, verbose=FALSE)
    skip_on_bioc(vdiffr::expect_doppelganger("Density Main Plot", grid::grid.draw(buildDensityPlot.out)))
})

test_that("buildDensityPlot is able to add layers to the plot", {
    
    skip_on_travis()
    
    test_plotBLayers <- list(ggplot2::geom_hline(yintercept=c(2e-9), colour="black", size=2), ggplot2::geom_vline(xintercept=c(1.5e8), colour="black", size=2))
    buildDensityPlot.out <- buildDensityPlot(RainfallPrimaryData.out, plotBLayers=test_plotBLayers , verbose=FALSE)
    skip_on_bioc(vdiffr::expect_doppelganger("Density Plot add layer", grid::grid.draw(buildDensityPlot.out)))
})


test_that("buildDensityPlot warns if plotBLayers is not passed as a list", {
    
    test_plotBLayers <- ggplot2::geom_hline(yintercept=c(30), colour="black", size=2)
    expect_error(buildDensityPlot(RainfallPrimaryData.out, plotBLayers=test_plotBLayers , verbose=FALSE))
})

test_that("buildDensityPlot warns if plotBLayers does not contain valid ggplot2 layers", {
    
    test_plotBLayers <- list(c("THIS IS A TEST"))
    expect_warning(buildDensityPlot(RainfallPrimaryData.out, plotBLayers=test_plotBLayers , verbose=FALSE))
})

test_that("buildDensityPlot works in verbose mode", {
    expect_message(buildDensityPlot(RainfallPrimaryData.out, plotBLayers=NULL, verbose=TRUE))
})

###################### test RainfallPlots class contruction ####################

context("RainfallPlots Contructor")

RainfallPlots.out <- RainfallPlots(RainfallPrimaryData.out, palette=NULL, pointSize=NULL, plotALayers=NULL,
                                   plotBLayers=NULL, verbose=FALSE)

test_that("RainfallPlots constructor outputs a S4 class object", {
    
    expect_s4_class(RainfallPlots.out, "RainfallPlots")
})

################################################################################
################### test Rainfall constructor and associate functions ##########

context("Rainfall Final Plot")

################### arrangeRainfallPlot ########################################

test_that("arrangeRainfallPlot plots a base plot", {
    
    skip_on_travis()
    
    arrangeRainfallPlot.out <- arrangeRainfallPlot(RainfallPlots.out, sectionHeights=NULL, verbose=FALSE)
    skip_on_bioc(vdiffr::expect_doppelganger("Final Rainfall Base", grid::grid.draw(arrangeRainfallPlot.out)))
})

test_that("arrangeRainfallPlot can alter section heights", {
    
    skip_on_travis()
    
    arrangeRainfallPlot.out <- arrangeRainfallPlot(RainfallPlots.out, sectionHeights=c(1, 1), verbose=FALSE)
    skip_on_bioc(vdiffr::expect_doppelganger("Final Rainfall alter section hieghts", grid::grid.draw(arrangeRainfallPlot.out)))
})

test_that("arrangeRainfallPlot warns if section heights does not match the number of plot elements", {
    
    expect_warning(arrangeRainfallPlot(RainfallPlots.out, sectionHeights=c(1, 1, 1), verbose=FALSE))
})

test_that("arrangeRainfallPlot warns if section heights are not numeric", {
    
    expect_warning(arrangeRainfallPlot(RainfallPlots.out, sectionHeights=as.character(c(1, 1, 1)), verbose=FALSE))
})

test_that("arrangeRainfallPlot works in verbose mode", {
    
    expect_message(arrangeRainfallPlot(RainfallPlots.out, sectionHeights=NULL, verbose=TRUE))
})

################################################################################
################## test Rainfall constructor and accessors #####################

Rainfall.out <- suppressWarnings(Rainfall(vepObject, BSgenome=BSgenome, palette=NULL, sectionHeights=NULL, chromosomes=c("chr1", "chr2"),
                                          sample=NULL, pointSize=NULL, verbose=FALSE, plotALayers=NULL, plotBLayers=NULL))

test_that("Rainfall constructor outputs a S4 class object", {
    
    expect_s4_class(Rainfall.out, "Rainfall")
})

############################ accessors #########################################

test_that("drawPlot constructs a Rainfall plot from grob objects in the Rainfall object", {
    
    skip_on_travis()
    
    skip_on_bioc(vdiffr::expect_doppelganger("drawPlot Rainfall", drawPlot(Rainfall.out)))
})

test_that("getData outputs error if no name or index is given", {
    
    expect_error(getData(Rainfall.out))
    
})

test_that("getData outputs error if index exceeds the number of slots", {
    
    expect_error(getData(Rainfall.out, index=10))
    
})

test_that("getData outputs error if supplied name is not a valid slot name", {
    
    expect_error(getData(Rainfall.out, name="shouldNotexist"))
    
})

test_that("getData retrieves specified slot data correctly", {
    
    expect_s3_class(getData(Rainfall.out, index=1), "data.table")
    expect_equivalent(getData(Rainfall.out, name="primaryData"), getData(Rainfall.out, index=1))
})

