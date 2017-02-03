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

