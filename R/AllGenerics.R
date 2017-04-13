#' Method getPosition
#' 
#' @name getPosition
#' @rdname getPosition-methods
#' @exportMethod getPosition
setGeneric(
    name="getPosition",
    def=function(object, ...){standardGeneric("getPosition")}
)

#' Method getMutation
#' 
#' @name getMutation
#' @rdname getMutation-methods
#' @exportMethod getMutation
setGeneric(
    name="getMutation",
    def=function(object, ...){standardGeneric("getMutation")}
)

#' Method getSample
#' 
#' @name getSample
#' @rdname getSample-methods
#' @exportMethod getSample
setGeneric(
    name="getSample",
    def=function(object, ...){standardGeneric("getSample")}
)

#' Method getMeta
#' 
#' @name getMeta
#' @rdname getMeta-methods
#' @exportMethod getMeta
setGeneric(
    name="getMeta",
    def=function(object, ...){standardGeneric("getMeta")}
)

#' Method toWaterfall
#'
#' @name toWaterfall
#' @rdname toWaterfall-methods
#' @noRd
setGeneric(
    name="toWaterfall",
    def=function(object, labelColumn, verbose, ...){standardGeneric("toWaterfall")}
)

#' Method sampSubset
#'
#' @name sampSubset
#' @rdname sampSubset-methods
#' @noRd
setGeneric(
    name="sampSubset",
    def=function(object, samples, verbose, ...){standardGeneric("sampSubset")}
)

#' Method calcSimpleMutationBurden
#'
#' @name calcSimpleMutationBurden
#' @rdname calcSimpleMutationBurden-methods
#' @noRd
setGeneric(
    name="calcSimpleMutationBurden",
    def=function(object, coverage, verbose, ...){standardGeneric("calcSimpleMutationBurden")}
)

#' Method calcComplexMutationBurden
#'
#' @name calcComplexMutationBurden
#' @rdname calcComplexMutationBurden-methods
#' @noRd
setGeneric(
    name="calcComplexMutationBurden",
    def=function(object, coverage, verbose, ...){standardGeneric("calcComplexMutationBurden")}
)

#' Method rmvSilentMutation
#'
#' @name rmvSilentMutation
#' @rdname rmvSilentMutation-methods
#' @noRd
setGeneric(
    name="rmvSilentMutation",
    def=function(object, verbose, ...){standardGeneric("rmvSilentMutation")}
)

#' Method geneSubset
#'
#' @name geneSubset
#' @rdname geneSubset-methods
#' @noRd
setGeneric(
    name="geneSubset",
    def=function(object, genes, verbose, ...){standardGeneric("geneSubset")}
)

#' Method mutHierarchySubset
#'
#' @name mutHierarchySubset
#' @rdname mutHierarchySubset-methods
#' @noRd
setGeneric(
    name="mutHierarchySubset",
    def=function(object, verbose, ...){standardGeneric("mutHierarchySubset")}
)

#' Method setMutationHierarchy
#'
#' @name setMutationHierarchy
#' @rdname setMutationHierarchy-methods
#' @noRd
setGeneric(
    name="setMutationHierarchy",
    def=function(object, mutationHierarchy, verbose, ...){standardGeneric("setMutationHierarchy")}
)

#' Method recurrenceSubset
#'
#' @name recurrenceSubset
#' @rdname recurrenceSubset-methods
#' @noRd
setGeneric(
    name="recurrenceSubset",
    def=function(object, recurrence, verbose, ...){standardGeneric("recurrenceSubset")}
)

#' Method orderGenes
#'
#' @name orderGenes
#' @rdname orderGenes-methods
#' @noRd
setGeneric(
    name="orderGenes",
    def=function(object, geneOrder, verbose, ...){standardGeneric("orderGenes")}
)

#' Method orderGenes
#'
#' @name maxGeneSubset
#' @rdname maxGeneSubset-methods
#' @noRd
setGeneric(
    name="maxGeneSubset",
    def=function(object, geneMax, verbose, ...){standardGeneric("maxGeneSubset")}
)

#' Method orderSamples
#'
#' @name orderSamples
#' @rdname orderSamples-methods
#' @noRd
setGeneric(
    name="orderSamples",
    def=function(object, sampleOrder, verbose, ...){standardGeneric("orderSamples")}
)

#' Method buildMutationPlot
#'
#' @name buildMutationPlot
#' @rdname buildMutationPlot-methods
#' @noRd
setGeneric(
    name="buildMutationPlot",
    def=function(object, plotA, plotATally, plotALayers, verbose, ...){standardGeneric("buildMutationPlot")}
)

#' Method constructGeneData
#'
#' @name constructGeneData
#' @rdname constructGeneData-methods
#' @noRd
setGeneric(
    name="constructGeneData",
    def=function(object, verbose, ...){standardGeneric("constructGeneData")}
)

#' Method buildGenePlot
#'
#' @name buildGenePlot
#' @rdname buildGenePlot-methods
#' @noRd
setGeneric(
    name="buildGenePlot",
    def=function(object, plotB, plotBTally, plotBLayers, verbose, ...){standardGeneric("buildGenePlot")}
)

#' Method buildWaterfallPlot
#'
#' @name buildWaterfallPlot
#' @rdname buildWaterfallPlot-methods
#' @noRd
setGeneric(
    name="buildWaterfallPlot",
    def=function(object, gridOverlay, drop, labelSize, labelAngle, sampleNames, xTitle, verbose, ...){standardGeneric("buildWaterfallPlot")}
)

#' Method formatClinicalData
#'
#' @name formatClinicalData
#' @rdname formatClinicalData-methods
#' @noRd
setGeneric(
    name="formatClinicalData",
    def=function(object, verbose, ...){standardGeneric("formatClinicalData")}
)

#' Method setClinicalPlotLayers
#'
#' @name setClinicalPlotLayers
#' @rdname setClinicalPlotLayers-methods
#' @noRd
setGeneric(
    name="setClinicalPlotLayers",
    def=function(object, legendColumns, palette, clinicalLayers, verbose, ...){standardGeneric("setClinicalPlotLayers")}
)

#' Method buildClinicalPlot
#'
#' @name buildClinicalPlot
#' @rdname buildClinicalPlot-methods
#' @noRd
setGeneric(
    name="buildClinicalPlot",
    def=function(object, verbose, ...){standardGeneric("buildClinicalPlot")}
)

#' Method getData
#' 
#' @name getData
#' @rdname getData-methods
#' @exportMethod getData
setGeneric(
    name="getData",
    def=function(object, ...){standardGeneric("getData")}
)

#' Method arrangeWaterfallPlot
#' 
#' @name arrangeWaterfallPlot
#' @rdname arrangeWaterfallPlot-methods
#' @exportMethod arrangeWaterfallPlot
setGeneric(
    name="arrangeWaterfallPlot",
    def=function(object, ...){standardGeneric("arrangeWaterfallPlot")}
)