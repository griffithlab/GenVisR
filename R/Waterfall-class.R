#' Class Waterfall
#' 
#' An S4 class for the waterfall plot object
#' @name Waterfall-class
#' @rdname Waterfall-class
#' @slot plotA gtable object for the top sub-plot.
#' @slot plotB gtable object for the left sub-plot.
#' @slot plotC gtable object for the main plot.
#' @slot plotD gtable object for the bottom sub-plot.
#' @slot Grob gtable object for the arranged plot.
#' @slot primaryData data.table object storing the primary data, should have
#' column names sample, gene, mutation, label.
#' @slot simpleMutationCounts data.table object storing simplified mutation
#' counts, should have column names sample, mutation, Freq, mutationBurden
#' @slot complexMutationData data.table object storing mutation counts per
#' mutation type should have column names sample, mutation, Freq, mutationBurden.
#' @slot GeneData
#' @slot ClinicalData
#' @exportClass Waterfall
#' @import methods
#' @importFrom gtable gtable
#' @importFrom data.table data.table

methods::setOldClass("gtable")
setClass("Waterfall",
         representation=representation(plotA="gtable",
                                       plotB="gtable",
                                       plotC="gtable",
                                       plotD="gtable",
                                       Grob="gtable",
                                       primaryData="data.table",
                                       simpleMutationCounts="data.table",
                                       complexMutationCounts="data.table",
                                       GeneData="data.table",
                                       ClinicalData="data.table"),
         validity=function(object){
             cat("!!!!! Waterfall~Inspector !!!!!\n")
         }
         )

#' Initalizer method for the Waterfall class
#' 
#' @name Waterfall
#' @rdname Waterfall-class
setMethod(f="initialize",
          signature="Waterfall",
          definition=function(.Object, input, labelColumn, samples, coverage,
                              noSynonymous, genes, verbose){

              # convert to waterfall format
              .Object@primaryData <- toWaterfall(input, labelColumn, verbose)
              
              # subset samples if specified
              .Object@primaryData <- sampSubset(.Object, samples, verbose)
              
              # calculate the mutation burden
              .Object@simpleMutationCounts <- calcSimpleMutationBurden(.Object, coverage, verbose)
              .Object@complexMutationCounts <- calcComplexMutationBurden(.Object, coverage, verbose)
              
              # remove silent mutations if specified
              if(noSynonymous) .Object@primaryData <- rmvSilentMutation(.Object, verbose)
              
              # subset on genes if specified
              .Object@primaryData <- geneSubset(.Object, genes, verbose)
              
              browser()
              return(.Object)
          })

#' Constructor for the Waterfall class.
#' 
#' @name Waterfall
#' @rdname Waterfall-class
#' @param input MutationAnnotationFormat class holding genomic information.
#' @param samples Character vector specifying samples to plot. If not NULL
#' all samples in "input" not specified with this parameter are removed.
#' @param coverage Integer specifying the size in base pairs of the genome
#' covered by sequence data from which mutations could be called. Required for
#' the mutation burden sub-plot (see details and vignette).
#' @param noSynonymous Boolean specifying if silent mutations should be removed
#' from the plot, this will not affect the mutation burden calculation.
#' @param verbose Boolean specifying if status messages should be reported
#' @export
Waterfall <- function(input, labelColumn=NULL, samples=NULL, coverage=NULL,
                      noSynonymous=FALSE, genes=NULL, verbose=FALSE){
    cat("!!!!! Waterfall~Constructor !!!!!\n")
    new("Waterfall", input=input, labelColumn=labelColumn, samples=samples, coverage=coverage,
        noSynonymous=noSynonymous, genes=genes, verbose=verbose)
}

#' @rdname Waterfall-methods
#' @aliases sampSubset,Waterfall
#' @param object Object of class waterfall
#' @param samples Character vector of samples to keep
#' @param verbose Boolean for status updates
#' @return data.table object subset on "samples" if "samples" is not null.
#' @noRd
setMethod(f="sampSubset",
          signature="Waterfall",
          definition=function(object, samples, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData
              
              # Dont do anything if samples is null
              if(is.null(samples)) return(primaryData)
              
              # print status message
              if(verbose){
                  memo <- paste("Restricting input to specified samples")
                  message(memo)
              }
              
              # make sure the samples variable is as expected
              if(class(samples) != 'character')
              {
                  memo <- paste0("argument supplied to samples is not a ",
                                 "character vector, attempting to coerce")
                  warning(memo)
                  samples <- as.character(samples)
              }
              
              # add in samples not in the original data frame x
              if(!all(samples %in% primaryData$sample))
              {
                  newSamples <- samples[!samples %in% primaryData$sample]
                  memo <- paste("The following samples were specified but not",
                                "found:", toString(newSamples), "adding these",
                               "samples to the input data!")
                  warning(memo)
                  primaryData <- rbind(primaryData,
                                       data.table("sample"=newSamples), fill=TRUE)
              }
              
              # remove any samples not requested to be plotted
              primaryData <- primaryData[primaryData$sample %in% samples,]
              primaryData$sample <- factor(primaryData$sample, levels=samples)
              
              return(primaryData)
          })

#' @rdname Waterfall-methods
#' @aliases calcSimpleMutationBurden,Waterfall
#' @param object Object of class waterfall
#' @param coverage Integer specifying the size in base pairs of genome from
#' which mutations could have been called (denominator of mutation burden).
#' @param verbose Boolean for status updates.
#' @return data.table object with simplified mutation counts
#' @noRd
#' @importFrom data.table setDT
#' @importFrom data.table data.table
setMethod(f="calcSimpleMutationBurden",
          signature="Waterfall",
          definition=function(object, coverage, verbose, ...){
              # if coverage is not specified to not return a mutation burden calculation
              if(!is.numeric(coverage)) return(data.table::data.table())
              
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData
              
              # status message
              if(verbose){
                  memo <- paste("Calculating a simplified mutation burden.")
                  message(memo)
              }
              
              # change mutation calls to either synonymous or non synonymous
              primaryData$simpleMutation <- NA
              silentMutations <- c("synonymous_variant", "silent", "synonymous")
              primaryData$simpleMutation[!toupper(primaryData$mutation) %in% toupper(silentMutations)] <- 'Non Synonymous'
              primaryData$simpleMutation[toupper(primaryData$mutation) %in% toupper(silentMutations)] <- 'Synonymous'
              primaryData$simpleMutation[is.na(primaryData$mutation)] <- NA
              primaryData$simpleMutation <- factor(primaryData$simpleMutation, levels=c('Synonymous', 'Non Synonymous'))
              
              # obtain a data table of mutation counts on the sample level
              simpleMutationCounts <- as.data.frame(table(primaryData[,c('sample', 'simpleMutation')]))
              data.table::setDT(simpleMutationCounts)
              colnames(simpleMutationCounts) <- c("sample", "mutation", "Freq")
              
              # mutation burden calculation
              simpleMutationCounts$mutationBurden <- simpleMutationCounts$Freq/coverage * 1000000
              return(simpleMutationCounts)
          })

#' @rdname Waterfall-methods
#' @aliases calcComplexMutationBurden,Waterfall
#' @param object Object of class waterfall
#' @param coverage Integer specifying the size in base pairs of genome from
#' which mutations could have been called (denominator of mutation burden).
#' @param verbose Boolean for status updates.
#' @return data.table object with complex mutation counts
#' @noRd
#' @importFrom data.table setDT
#' @importFrom data.table data.table
setMethod(f="calcComplexMutationBurden",
          signature="Waterfall",
          definition=function(object, coverage, verbose, ...){
              # if coverage is not specified correctly do not return a mutation burden calculation
              if(!is.numeric(coverage)) return(data.table::data.table())
              
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData
              
              # status message
              if(verbose){
                  memo <- paste("Calculating a complex mutation burden.")
                  message(memo)
              }
              
              # obtain a data table of mutation counts on the sample level
              complexMutationCounts <- as.data.frame(table(primaryData[,c('sample', 'mutation')]))
              data.table::setDT(complexMutationCounts)
              
              # mutation burden calculation
              complexMutationCounts$mutationBurden <- complexMutationCounts$Freq/coverage * 1000000
              return(complexMutationCounts)
          })

#' @rdname Waterfall-methods
#' @aliases rmvSilentMutation,Waterfall
#' @param object Object of class waterfall
#' @param verbose Boolean for status updates
#' @return data.table object with silent mutations removed from primaryData slot
#' @noRd
setMethod(f="rmvSilentMutation",
          signature="Waterfall",
          definition=function(object, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData
              
              # store how many mutations there are originaly to figure out how
              # silent mutations are removed
              totalMut <- length(na.omit(primaryData$mutation))
              
              # define what a silent mutation is
              silentMutations <- c("synonymous", "silent", "synonymous_variant")
              
              # assign NA to any rows/columns with silent mutations, keep the
              # sample so the cohort size is not accidently reduced
              primaryData[toupper(primaryData$trv_type) == toupper(silentMutations), c('gene', 'mutation', 'label')] <- NA
              
              # figure out number of silent mutations that existed
              silentMutCount <- totalMut - length(na.omit(primaryData$mutation))
              
              # print status message
              if(verbose){
                  memo <- paste("Found", silentMutCount,"silent mutations...",
                                "removing")
                  message(memo)
              }
              
              return(primaryData)
          })

#' @rdname Waterfall-methods
#' @aliases geneSubset,Waterfall
#' @param object Object of class waterfall
#' @param verbose Boolean for status updates
#' @return data.table object subset on gene if gene is not NULL
#' @noRd
setMethod(f="geneSubset",
          signature="Waterfall",
          definition=function(object, genes, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData
              
              # Dont do anything if samples is null
              if(is.null(samples)) return(primaryData)
              
              # print status message
              if(verbose){
                  memo <- paste("Restricting input to specified genes")
                  message(memo)
              }
              
              # Perform quality checks
              if(class(genes) != 'character')
              {
                  memo <- paste("argument supplied to genes is not of class",
                                 "character, attempting to coerce!")
                  warning(memo)
                  genes <- as.character(genes)
              }
              
              # check if a gene was specified but is not in the data
              if(!all(genes %in% primaryData$gene))
              {
                  newGenes <- genes[!genes %in% primaryData$gene]
                  memo <- paste("The following genes were specified but were",
                                "not found in the data:", toString(newGenes),"
                                this may be because they were filtered in a", 
                                "previous step... adding these to the input.")
                  warning(memo)
                  primaryData <- rbind(primaryData,
                                       data.table("gene"=newGenes), fill=TRUE)
              }
              
              # keep all genes specified and any NA so samples are kept
              genes <- c(genes, NA)
              primaryData <- primaryData[primaryData$gene %in% genes,]
              
              return(primaryData)
          })
