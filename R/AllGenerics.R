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
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="annoGenomeCoord",
    def=function(object, ...){standardGeneric("annoGenomeCoord")}
)

#' Method chrSubset
#'
#' @name chrSubset
#' @rdname chrSubset-methods
#' @param ... additional arguments to passed
#' @noRd
setGeneric(
    name="chrSubset",
    def=function(object, ...){standardGeneric("chrSubset")}
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