#' Method getPosition
#' 
#' @name getPosition
#' @rdname getPosition-methods
#' @param object Object of class VEP, GMS, or MutationAnnotationFormat
#' @param ... additional arguments to passed
#' @exportMethod getPosition
setGeneric(
    name="getPosition",
    def=function(object, ...){standardGeneric("getPosition")}
)

#' Method getMutation
#' 
#' @name getMutation
#' @rdname getMutation-methods
#' @param object Object of class VEP, GMS, or MutationAnnotationFormat
#' @param ... additional arguments to passed
#' @exportMethod getMutation
setGeneric(
    name="getMutation",
    def=function(object, ...){standardGeneric("getMutation")}
)

#' Method getSample
#' 
#' @name getSample
#' @rdname getSample-methods
#' @param object Object of class VEP, GMS, or MutationAnnotationFormat
#' @param ... additional arguments to passed
#' @exportMethod getSample
setGeneric(
    name="getSample",
    def=function(object, ...){standardGeneric("getSample")}
)

#' Method getMeta
#' 
#' @name getMeta
#' @rdname getMeta-methods
#' @param object Object of class VEP, GMS, or MutationAnnotationFormat
#' @param ... additional arguments to passed
#' @exportMethod getMeta
setGeneric(
    name="getMeta",
    def=function(object, ...){standardGeneric("getMeta")}
)

#' Method toWaterfall
#'
#' @name toWaterfall
#' @rdname toWaterfall-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="toWaterfall",
    def=function(object, labelColumn, verbose, ...){standardGeneric("toWaterfall")}
)

#' Method sampSubset
#'
#' @name sampSubset
#' @rdname sampSubset-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="sampSubset",
    def=function(object, samples, verbose, ...){standardGeneric("sampSubset")}
)

#' Method calcSimpleMutationBurden
#'
#' @name calcSimpleMutationBurden
#' @rdname calcSimpleMutationBurden-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="calcSimpleMutationBurden",
    def=function(object, coverage, verbose, ...){standardGeneric("calcSimpleMutationBurden")}
)

#' Method calcComplexMutationBurden
#'
#' @name calcComplexMutationBurden
#' @rdname calcComplexMutationBurden-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="calcComplexMutationBurden",
    def=function(object, coverage, verbose, ...){standardGeneric("calcComplexMutationBurden")}
)

#' Method rmvMutation
#'
#' @name rmvMutation
#' @rdname rmvMutation-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="rmvMutation",
    def=function(object, verbose, ...){standardGeneric("rmvMutation")}
)

#' Method geneSubset
#'
#' @name geneSubset
#' @rdname geneSubset-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="geneSubset",
    def=function(object, genes, verbose, ...){standardGeneric("geneSubset")}
)

#' Method mutHierarchySubset
#'
#' @name mutHierarchySubset
#' @rdname mutHierarchySubset-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="mutHierarchySubset",
    def=function(object, verbose, ...){standardGeneric("mutHierarchySubset")}
)

#' Method setMutationHierarchy
#'
#' @name setMutationHierarchy
#' @rdname setMutationHierarchy-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="setMutationHierarchy",
    def=function(object, mutationHierarchy, verbose, ...){standardGeneric("setMutationHierarchy")}
)

#' Method recurrenceSubset
#'
#' @name recurrenceSubset
#' @rdname recurrenceSubset-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="recurrenceSubset",
    def=function(object, recurrence, verbose, ...){standardGeneric("recurrenceSubset")}
)

#' Method orderGenes
#'
#' @name orderGenes
#' @rdname orderGenes-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="orderGenes",
    def=function(object, geneOrder, verbose, ...){standardGeneric("orderGenes")}
)

#' Method orderGenes
#'
#' @name maxGeneSubset
#' @rdname maxGeneSubset-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="maxGeneSubset",
    def=function(object, geneMax, verbose, ...){standardGeneric("maxGeneSubset")}
)

#' Method orderSamples
#'
#' @name orderSamples
#' @rdname orderSamples-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="orderSamples",
    def=function(object, sampleOrder, verbose, ...){standardGeneric("orderSamples")}
)

#' Method buildMutationPlot
#'
#' @name buildMutationPlot
#' @rdname buildMutationPlot-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="buildMutationPlot",
    def=function(object, plotA, plotATally, plotALayers, verbose, ...){standardGeneric("buildMutationPlot")}
)

#' Method constructGeneData
#'
#' @name constructGeneData
#' @rdname constructGeneData-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="constructGeneData",
    def=function(object, verbose, ...){standardGeneric("constructGeneData")}
)

#' Method buildGenePlot
#'
#' @name buildGenePlot
#' @rdname buildGenePlot-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="buildGenePlot",
    def=function(object, plotB, plotBTally, plotBLayers, verbose, ...){standardGeneric("buildGenePlot")}
)

#' Method buildWaterfallPlot
#'
#' @name buildWaterfallPlot
#' @rdname buildWaterfallPlot-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="buildWaterfallPlot",
    def=function(object, gridOverlay, drop, labelSize, labelAngle, sampleNames, xTitle, verbose, ...){standardGeneric("buildWaterfallPlot")}
)

#' Method formatClinicalData
#'
#' @name formatClinicalData
#' @rdname formatClinicalData-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="formatClinicalData",
    def=function(object, verbose, ...){standardGeneric("formatClinicalData")}
)

#' Method setClinicalPlotLayers
#'
#' @name setClinicalPlotLayers
#' @rdname setClinicalPlotLayers-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="setClinicalPlotLayers",
    def=function(object, legendColumns, palette, clinicalLayers, verbose, ...){standardGeneric("setClinicalPlotLayers")}
)

#' Method buildClinicalPlot
#'
#' @name buildClinicalPlot
#' @rdname buildClinicalPlot-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="buildClinicalPlot",
    def=function(object, verbose, ...){standardGeneric("buildClinicalPlot")}
)

#' Method getData
#' 
#' @name getData
#' @rdname getData-methods
#' @param object Object of class Clinical, 
#' @param name String corresponding to the slot for which to extract data from.
#' @param index Integer specifying the slot for which to extract data from.
#' @param ... additional arguments to passed
#' @details The getData method is an accessor function used to access data held
#' in GenVisR objects.
#' @exportMethod getData
setGeneric(
    name="getData",
    def=function(object, ...){standardGeneric("getData")}
)

#' Method getVersion
#' 
#' @name getVersion
#' @rdname getVersion-methods
#' @param object Object of class VEP, GMS, or MutationAnnotationFormat
#' @param ... additional arguments to passed
#' @exportMethod getVersion
setGeneric(
    name="getVersion",
    def=function(object, ...){standardGeneric("getVersion")}
)

#' Method getPath
#' 
#' @name getPath
#' @rdname getPath-methods
#' @param object Object of class VEP, GMS, or MutationAnnotationFormat
#' @param ... additional arguments to passed
#' @exportMethod getPath
setGeneric(
    name="getPath",
    def=function(object, ...){standardGeneric("getPath")}
)

#' Method getLayers
#' 
#' @name getLayers
#' @rdname getLayers-methods
#' @param object Object of class Clinical
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="getLayers",
    def=function(object, ...){standardGeneric("getLayers")}
)

#' Method arrangeWaterfallPlot
#' 
#' @name arrangeWaterfallPlot
#' @rdname arrangeWaterfallPlot-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="arrangeWaterfallPlot",
    def=function(object, ...){standardGeneric("arrangeWaterfallPlot")}
)

#' Method drawPlot
#' 
#' @name drawPlot
#' @rdname drawPlot-methods
#' @param object Object of class Waterfall, MutSpectra, or Clinical
#' @param ... additional arguments to passed
#' @details The drawPlot method is used to draw plots created by GenVisR plot 
#' constructor functions.
#' @exportMethod drawPlot
setGeneric(
    name="drawPlot",
    def=function(object, ...){standardGeneric("drawPlot")}
)

#' Method parseHeader
#' 
#' @name parseHeader
#' @rdname parseHeader-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="parseHeader",
    def=function(object, ...){standardGeneric("parseHeader")}
)

#' Method parseDescription
#' 
#' @name parseDescription
#' @rdname parseDescription-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="parseDescription",
    def=function(object, ...){standardGeneric("parseDescription")}
)

#' Method parseExtra
#' 
#' @name parseExtra
#' @rdname parseExtra-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="parseExtra",
    def=function(object, ...){standardGeneric("parseExtra")}
)

#' Method writeData
#' 
#' @name writeData
#' @rdname writeData-methods
#' @param object Object of class VEP
#' @param file Character string specifying a file to send output to.
#' @param sep Delimiter used when writing output, defaults to tab.
#' @param ... additional arguments to passed
#' @details The writeData method is used to output data held in GenVisR objects
#' to a file.
#' @exportMethod writeData
setGeneric(
    name="writeData",
    def=function(object, ...){standardGeneric("writeData")}
)

#' Method geneFilter
#' 
#' @name geneFilter
#' @rdname geneFilter-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="geneFilter",
    def=function(object, ...){standardGeneric("geneFilter")}
)

#' Method toMutSpectra
#' 
#' @name toMutSpectra
#' @rdname toMutSpectra
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="toMutSpectra",
    def=function(object, ...){standardGeneric("toMutSpectra")}
)

#' Method annoMutSpectra
#' 
#' @name annoMutSpectra
#' @rdname annoMutSpectra
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="annoMutSpectra",
    def=function(object, ...){standardGeneric("annoMutSpectra")}
)

#' Method calcMutSpectra
#' 
#' @name calcMutSpectra
#' @rdname calcMutSpectra
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="calcMutSpectra",
    def=function(object, ...){standardGeneric("calcMutSpectra")}
)

#' Method sortSamples
#' 
#' @name sortSamples
#' @rdname sortSamples
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="sortSamples",
    def=function(object, ...){standardGeneric("sortSamples")}
)

#' Method buildFrequencyPlot
#' 
#' @name buildFrequencyPlot
#' @rdname buildFrequencyPlot
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="buildFrequencyPlot",
    def=function(object, ...){standardGeneric("buildFrequencyPlot")}
)

#' Method buildProportionPlot
#' 
#' @name buildProportionPlot
#' @rdname buildProportionPlot
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="buildProportionPlot",
    def=function(object, ...){standardGeneric("buildProportionPlot")}
)

#' Method arrangeMutSpectraPlot
#' 
#' @name arrangeMutSpectraPlot
#' @rdname arrangeMutSpectraPlot
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="arrangeMutSpectraPlot",
    def=function(object, ...){standardGeneric("arrangeMutSpectraPlot")}
)

#' Method getGrob
#' 
#' @name getGrob
#' @rdname getGrob-methods
#' @param object Object of clas MutSpectra
#' @param index integer specifying the plot index to extract
#' @param ... additional arguments to passed
#' @exportMethod getGrob
setGeneric(
    name="getGrob",
    def=function(object, ...){standardGeneric("getGrob")}
)

#' Method getDescription
#' 
#' @name getDescription
#' @rdname getDescription-methods
#' @param object Object of class VEP
#' @param ... additional arguments to passed
#' @exportMethod getDescription
setGeneric(
    name="getDescription",
    def=function(object, ...){standardGeneric("getDescription")}
)

#' Method getHeader
#' 
#' @name getHeader
#' @rdname getHeader-methods
#' @param object Object of class VEP
#' @param ... additional arguments to passed
#' @exportMethod getHeader
setGeneric(
    name="getHeader",
    def=function(object, ...){standardGeneric("getHeader")}
)

################################################################################
##### Functions used for lohSpec ###############################################
################################################################################
#' Method getVarScan
#' 
#' @name getVarScan
#' @rdname getVarScan-methods
#' @param object Object of class VarScanFormat
#' @param ... additional arguments to passed
#' @exportMethod getVarScan
setGeneric(
    name="getVarScan",
    def=function(object, ...){standardGeneric("getVarScan")}
)

#' Method getLohData
#' 
#' @name getLohData
#' @rdname getLohData-methods
#' @param object Object of class VarScanFormat
#' @param ... additional arguments to passed
#' @exportMethod getLohData
setGeneric(
    name="getLohData",
    def=function(object, ...){standardGeneric("getLohData")}
)

#' Method lohSpec_qual
#' 
#' @name lohSpec_qual
#' @rdname lohSpec_qual-methods
#' @param object Object of class VarScanFormat
#' @param ... additional arguments to passed
setGeneric(
    name="lohSpec_qual",
    def=function(object, ...){standardGeneric("lohSpec_qual")}
)

#' Method toRainfall
#'
#' @name toRainfall
#' @rdname toRainfall-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="toRainfall",
    def=function(object, ...){standardGeneric("toRainfall")}
)

#' Method annoRainfall
#'
#' @name annoRainfall
#' @rdname annoRainfall-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="annoRainfall",
    def=function(object, ...){standardGeneric("annoRainfall")}
)

#' Method calcMutDist
#'
#' @name calcMutDist
#' @rdname calcMutDist-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="calcMutDist",
    def=function(object, ...){standardGeneric("calcMutDist")}
)

#' Method annoGenomeCoord
#'
#' @name annoGenomeCoord
#' @rdname annoGenomeCoord-methods
#' @param object object of class data.table
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="annoGenomeCoord",
    def=function(object, ...){standardGeneric("annoGenomeCoord")}
)

#' Method annoGenomeCoordSv
#'
#' @name annoGenomeCoordSv
#' @rdname annoGenomeCoordSv-methods
#' @param object object of class data.table
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="annoGenomeCoordSv",
    def=function(object, ...){standardGeneric("annoGenomeCoordSv")}
)

#' Method getLohSlidingWindow
#' 
#' @name getLohSlidingWindow
#' @rdname getLohSlidingWindow-methods
#' @param object Object of class VarScanFormat
#' @param ... additional arguments to passed
setGeneric(
    name="getLohSlidingWindow",
    def=function(object, ...){standardGeneric("getLohSlidingWindow")}
)

#' Method getLohCalculation
#' 
#' @name getLohCalculation
#' @rdname getLohCalculation-methods
#' @param object Object of class VarScanFormat
#' @param ... additional arguments to passed
setGeneric(
    name="getLohCalculation",
    def=function(object, ...){standardGeneric("getLohCalculation")}
)

#' Method getLohStepCalculation
#' 
#' @name getLohStepCalculation
#' @rdname getLohStepCalculation-methods
#' @param object Object of class VarScanFormat
#' @param ... additional arguments to passed
setGeneric(
    name="getLohStepCalculation",
    def=function(object, ...){standardGeneric("getLohStepCalculation")}
)

#' Method getLohSegmentation
#' 
#' @name getLohSegmentation
#' @rdname getLohSegmentation-methods
#' @param object Object of class VarScanFormat
#' @param ... additional arguments to passed
setGeneric(
    name="getLohSegmentation",
    def=function(object, ...){standardGeneric("getLohSegmentation")}
)

#' Method getLohFreq
#' 
#' @name getLohFreq
#' @rdname getLohFreq-methods
#' @param object Object of class VarScanFormat
#' @param ... additional arguments to passed
setGeneric(
    name="getLohFreq",
    def=function(object, ...){standardGeneric("getLohFreq")}
)

#' Method buildLohFreq
#' 
#' @name buildLohFreq
#' @rdname buildLohFreq-methods
#' @param object Object of class VarScanFormat
#' @param ... additional arguments to passed
setGeneric(
    name="buildLohFreq",
    def=function(object, ...){standardGeneric("buildLohFreq")}
)

#' Method lohSpec_buildMainPlot
#' 
#' @name lohSpec_buildMainPlot
#' @rdname lohSpec_buildMainPlot-methods
#' @param object Object of class VarScanFormat
#' @param ... additional arguments to passed
setGeneric(
    name="lohSpec_buildMainPlot",
    def=function(object, ...){standardGeneric("lohSpec_buildMainPlot")}
)

#' Method arrangeLohPlots
#' 
#' @name arrangeLohPlots
#' @rdname arrangeLohPlots-methods
#' @param object Object of class VarScanFormat
#' @param ... additional arguments to passed
setGeneric(
    name="arrangeLohPlots",
    def=function(object, ...){standardGeneric("arrangeLohPlots")}
)

#' Method chrSubset
#'
#' @name chrSubset
#' @rdname chrSubset-methods
#' @param object object of class data.table
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="chrSubset",
    def=function(object, ...){standardGeneric("chrSubset")}
)

#' Method chrSubsetSv
#'
#' @name chrSubsetSv
#' @rdname chrSubsetSv-methods
#' @param object object of class data.table
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="chrSubsetSv",
    def=function(object, ...){standardGeneric("chrSubsetSv")}
)

#' Method sampleSubset
#'
#' @name sampleSubset
#' @rdname sampleSubset-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="sampleSubset",
    def=function(object, ...){standardGeneric("sampleSubset")}
)

#' Method highlightSampleData
#'
#' @name formatSample
#' @rdname formatSample-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="formatSample",
    def=function(object, ...){standardGeneric("formatSample")}
)

#' Method buildRainfallPlot
#'
#' @name buildRainfallPlot
#' @rdname buildRainfallPlot-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="buildRainfallPlot",
    def=function(object, ...){standardGeneric("buildRainfallPlot")}
)

#' Method buildDensityPlot
#'
#' @name buildDensityPlot
#' @rdname buildDensityPlot-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="buildDensityPlot",
    def=function(object, ...){standardGeneric("buildDensityPlot")}
)

#' Method arrangeRainfallPlot
#'
#' @name arrangeRainfallPlot
#' @rdname arrangeRainfallPlot-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="arrangeRainfallPlot",
    def=function(object, ...){standardGeneric("arrangeRainfallPlot")}
)

#' Method getCnvData
#' 
#' @name getCnvData
#' @rdname getCnvData-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="getCnvData",
    def=function(object, ...){standardGeneric("getCnvData")}
)

#' Method getCnSegmentation
#' 
#' @name getCnSegmentation
#' @rdname getCnSegmentation-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="getCnSegmentation",
    def=function(object, ...){standardGeneric("getCnSegmentation")}
)

#' Method buildCnPlot
#' 
#' @name buildCnPlot
#' @rdname buildCnPlot-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="buildCnPlot",
    def=function(object, ...){standardGeneric("buildCnPlot")}
)

#' Method buildSomaticLohPlot
#' 
#' @name buildSomaticLohPlot
#' @rdname buildSomaticLohPlot-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="buildSomaticLohPlot",
    def=function(object, ...){standardGeneric("buildSomaticLohPlot")}
)

#' Method buildGermlineLohPlot
#' 
#' @name buildGermlineLohPlot
#' @rdname buildGermlineLohPlot-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="buildGermlineLohPlot",
    def=function(object, ...){standardGeneric("buildGermlineLohPlot")}
)

#' Method arrangeCnLohPlots
#' 
#' @name arrangeCnLohPlots
#' @rdname arrangeCnLohPlots-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="arrangeCnLohPlots",
    def=function(object, ...){standardGeneric("arrangeCnLohPlots")}
)

#' Method removeGapsSegmentation
#' 
#' @name removeGapsSegmentation
#' @rdname removeGapsSegmentation-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="removeGapsSegmentation",
    def=function(object, ...){standardGeneric("removeGapsSegmentation")}
)

#' Method getVcf
#' 
#' @name getVcf
#' @rdname getVcf-methods
#' @param ... additional arguments to passed
#' @exportMethod getSample
setGeneric(
    name="getVcf",
    def=function(object, ...){standardGeneric("getVcf")}
)

#' Method filterStructuralVariant
#' 
#' @name filterStructuralVariant
#' @rdname filterStructuralVariant-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="filterStructuralVariant",
    def=function(object, ...){standardGeneric("filterStructuralVariant")}
)

#' Method annotateSV
#' 
#' @name annotateSV
#' @rdname annotateSV-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="annotateSV",
    def=function(object, ...){standardGeneric("annotateSV")}
)

#' Method countGenes
#' 
#' @name countGenes
#' @rdname countGenes-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="countGenes",
    def=function(object, ...){standardGeneric("countGenes")}
)

#' Method getStructuralVariantWindow
#' 
#' @name getStructuralVariantWindow
#' @rdname getStructuralVariantWindow-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="getStructuralVariantWindow",
    def=function(object, ...){standardGeneric("getStructuralVariantWindow")}
)

#' Method buildSvPlot
#' 
#' @name buildSvPlot
#' @rdname buildSvPlot-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="buildSvPlot",
    def=function(object, ...){standardGeneric("buildSvPlot")}
)

#' Method svCytobands
#' 
#' @name svCytobands
#' @rdname svCytobands-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="svCytobands",
    def=function(object, ...){standardGeneric("svCytobands")}
)

#' Method adjustCentromeres
#' 
#' @name adjustCentromeres
#' @rdname adjustCentromeres-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="adjustCentromeres",
    def=function(object, ...){standardGeneric("adjustCentromeres")}
)

#' Method getVcfData
#' 
#' @name getVcfData
#' @rdname getVcfData-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="getVcfData",
    def=function(object, ...){standardGeneric("getVcfData")}
)

#' Method checkSvInputParameters
#' 
#' @name checkSvInputParameters
#' @rdname checkSvInputParameters-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="checkSvInputParameters",
    def=function(object, ...){standardGeneric("checkSvInputParameters")}
)