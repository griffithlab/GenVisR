#' Class Waterfall
#' 
#' An S4 class for the waterfall plot object
#' @name Waterfall-class
#' @rdname Waterfall-class
#' @slot PlotA gtable object for the top sub-plot.
#' @slot PlotB gtable object for the left sub-plot.
#' @slot PlotC gtable object for the main plot.
#' @slot PlotD gtable object for the bottom sub-plot.
#' @slot Grob gtable object for the arranged plot.
#' @slot primaryData data.table object storing the primary data, should have
#' column names sample, gene, mutation, label.
#' @slot simpleMutationCounts data.table object storing simplified mutation
#' counts, should have column names sample, mutation, Freq, mutationBurden
#' @slot complexMutationCounts data.table object storing mutation counts per
#' mutation type should have column names sample, mutation, Freq, mutationBurden.
#' @slot geneData data.table object storing gene counts, should have column
#' names gene, mutation, count.
#' @slot ClinicalData data.table object stroring the data used to plot the
#' clinical sub-plot.
#' @slot MutationHierarchy data.table object storing the hierarchy of mutation
#' type in order of most to least important and the mapping of mutation type to
#' color. Should have column names mutation, color, and label.
#' @exportClass Waterfall
#' @import methods
#' @importFrom gtable gtable
#' @importFrom data.table data.table

methods::setOldClass("gtable")
setClass("Waterfall",
         representation=representation(PlotA="gtable",
                                       PlotB="gtable",
                                       PlotC="gtable",
                                       PlotD="gtable",
                                       Grob="gtable",
                                       primaryData="data.table",
                                       simpleMutationCounts="data.table",
                                       complexMutationCounts="data.table",
                                       geneData="data.table",
                                       ClinicalData="data.table",
                                       MutationHierarchy="data.table"),
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
                              noSynonymous, genes, mutationHierarchy, recurrence, 
                              geneOrder, geneMax, sampleOrder, plotA,
                              plotATally, plotALayers, plotB, plotBTally,
                              plotBLayers, gridOverlay, drop, labelSize, labelAngle,
                              sampleNames, clinical, sectionHeights, sectionWidths, 
                              verbose){

              # convert to initial data to waterfall format
              .Object@primaryData <- toWaterfall(input, labelColumn, verbose)
              
              # assign the mapping of mutations and colors
              .Object@MutationHierarchy <- setMutationHierarchy(input, mutationHierarchy, verbose)
              
              # subset samples if specified
              .Object@primaryData <- sampSubset(.Object, samples, verbose)
              
              # calculate the frequency and mutation burden
              .Object@simpleMutationCounts <- calcSimpleMutationBurden(.Object, coverage, verbose)
              .Object@complexMutationCounts <- calcComplexMutationBurden(.Object, coverage, verbose)
              
              # remove silent mutations if specified
              if(noSynonymous) .Object@primaryData <- rmvSilentMutation(.Object, verbose)
              
              # subset on genes if specified
              .Object@primaryData <- geneSubset(.Object, genes, verbose)
              
              # remove entries for the same gene/sample based on a hierarchy leaving one
              .Object@primaryData <- mutHierarchySubset(.Object, verbose)
              
              # subset on recurrence of mutations
              .Object@primaryData <- recurrenceSubset(.Object, recurrence, verbose)
              
              # set the order of genes for plotting
              .Object@primaryData <- orderGenes(.Object, geneOrder, verbose)
              
              # limit to a maximum number of genes
              .Object@primaryData <- maxGeneSubset(.Object, geneMax, verbose)
              
              # set the order of samples for plotting
              .Object@primaryData <- orderSamples(.Object, sampleOrder, verbose)
              browser()
              # create the top sub-plot
              .Object@PlotA <- buildMutationPlot(.Object, plotA, plotATally, 
                                                 plotALayers, verbose)
              # summarize gene level data
              .Object@geneData <- constructGeneData(.Object, verbose)
              
              # create left sub-plot
              .Object@PlotB <- buildGenePlot(.Object, plotB, plotBTally,
                                             plotBLayers, verbose)
              
              # add the clinical data
              if(is.null(clinical)){
                  .Object@ClinicalData <- data.table::data.table()
              } else {
                  # there be dragons here, fix!
                  .Object@ClinicalData <- getData(clinical)
                  .Object@ClinicalData <- formatClinicalData(.Object, verbose)
              }

              # add the clinical data plot
              .Object@PlotD <- buildClinicalPlot(.Object, clinicalLayers=clinical@clinicalLayers)
              
              # create the main plot
              xTitle <- TRUE #tmporary remove when clinical object is defined
              .Object@PlotC <- buildWaterfallPlot(.Object, gridOverlay, drop,
                                                  labelSize, labelAngle, xTitle,
                                                  sampleNames, verbose)
              
              # align all plots together
              .Object@Grob <- arrangeWaterfallPlot(.Object, sectionHeights=sectionHeights,
                                                   sectionWidths=sectionWidths,
                                                   verbose=verbose)
              
              return(.Object)
          })

#' Constructor for the Waterfall class.
#' 
#' @name Waterfall
#' @rdname Waterfall-class
#' @param input Object of class MutationAnnotationFormat, or GMS.
#' @param labelColumn Character vector specifying a column name from which to
#' extract labels for cells.
#' @param samples Character vector specifying samples to plot. If not NULL
#' all samples in "input" not specified with this parameter are removed.
#' @param coverage Integer specifying the size in base pairs of the genome
#' covered by sequence data from which mutations could be called. Required for
#' the mutation burden sub-plot (see details and vignette).
#' @param noSynonymous Boolean specifying if silent mutations should be removed
#' from the plot, this will not affect the mutation burden calculation.
#' @param mutationHierarchy Data.table object with rows specifying the order of
#' mutations from most to least deleterious and column names "mutation" and
#' "color". Used to change the default colors and/or to give priority to a
#' mutation for the same gene/sample (see details and vignette).
#' @param recurrence Numeric value between 0 and 1 specifying a
#' mutation recurrence cutoff. Genes which do not have mutations in the
#' proportion of samples defined are removed.
#' @param geneOrder Character vector specifying the order in which to plot
#' genes.
#' @param geneMax Integer specifying the maximum number of genes to be plotted.
#' Genes kept will be choosen based on the reccurence of mutations in samples.
#' Unless geneOrder is specified.
#' @param sampOrder Character vector specifying the order in which to plot
#' samples.
#' @param plotA String specifying the type of plot for the top sub-plot, one of
#' "burden", "frequency", or NULL for a mutation burden (requires coverage to be
#' specified), frequency of mutations, or no plot respectively.
#' @param plotALayers list of ggplot2 layers to be passed to the plot.
#' @param plotATally String specifying one of "simple" or "complex" for a
#' simplified or complex tally of mutations respectively.
#' @param plotB String specifying the type of plot for the left sub-plot, one of
#' "proportion", "frequency", or NULL for a plot of gene proportions frequencies
#' , or no plot respectively.
#' @param plotBTally String specifying one of "simple" or "complex" for a
#' simplified or complex tally of genes respectively.
#' @param plotBLayers list of ggplot2 layers to be passed to the plot.
#' @param gridOverlay Boolean specifying if a grid should be overlayed on the
#' waterfall plot.
#' @param drop Boolean specifying if unused mutations should be dropped from the
#' legend.
#' @param labelSize Integer specifying the size of label text
#' @param labelAngle Numeric value specifying the angle of label text
#' @param sampleNames Boolean specifying if samples should be labeled on the plot.
#' @param clinical Object of class Clinical, used for adding a clinical data subplot.
#' @param sectionHeights Numeric vector specifying relative heights of each plot section,
#' should sum to one. Expects a value for each section.
#' @param sectionWidths Numeric vector specifying relative heights of each plot section,
#' should sum to one. Expects a value for each section.
#' @param verbose Boolean specifying if status messages should be reported
#' @export
Waterfall <- function(input, labelColumn=NULL, samples=NULL, coverage=NULL,
                      noSynonymous=FALSE, genes=NULL, mutationHierarchy=NULL,
                      recurrence=NULL, geneOrder=NULL, geneMax=NULL,
                      sampleOrder=NULL, plotA=c("frequency", "burden", NULL),
                      plotATally=c("simple", "complex"), plotALayers=NULL,
                      plotB=c("proportion", "frequency", NULL),
                      plotBTally=c("simple", "complex"), plotBLayers=NULL,
                      gridOverlay=FALSE, drop=TRUE, labelSize=5, labelAngle=0,
                      sampleNames=TRUE, clinical=NULL, sectionHeights=NULL,
                      sectionWidths=NULL, verbose=FALSE){
    cat("!!!!! Waterfall~Constructor !!!!!\n")
    new("Waterfall", input=input, labelColumn=labelColumn, samples=samples, coverage=coverage,
        noSynonymous=noSynonymous, genes=genes, mutationHierarchy=mutationHierarchy,
        recurrence=recurrence, geneOrder=geneOrder, geneMax=geneMax, sampleOrder=sampleOrder,
        plotA=plotA, plotATally=plotATally, plotALayers=plotALayers, plotB=plotB,
        plotBTally=plotBTally, plotBLayers=plotBLayers, gridOverlay=gridOverlay, drop=drop,
        labelSize=labelSize, labelAngle=labelAngle, sampleNames=sampleNames,
        clinical=clinical, sectionHeights=sectionHeights, sectionWidths=sectionWidths,
        verbose=verbose)
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
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData
              
              # quality checks
              if(length(coverage) > 1 && !is.null(coverage)){
                  memo <- paste("coverage has a length > 1, using only the",
                                "first element.")
                  warning(memo)
                  coverage <- coverage[1]
              }
              if(!is.numeric(coverage) && !is.null(coverage)){
                  memo <- paste("coverage is not numeric, attempting to coerce.")
                  warning(memo)
                  coverage <- as.numeric(coverage)
              }
              
              # status message
              if(verbose){
                  memo <- paste("Calculating frequency and mutation burden.")
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
              
              # if coverage is not specified return just frequencies
              if(!is.numeric(coverage)){
                  if(verbose){
                      memo <- paste("coverage not specified, could not",
                                    "calculate the mutation burden")
                      message(memo)
                      simpleMutationCounts$mutationBurden <- NA
                      return(simpleMutationCounts)
                  }
              } 
              
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
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData
              
              # quality checks
              if(length(coverage) > 1 && !is.null(coverage)){
                  memo <- paste("coverage has a length > 1, using only the",
                                "first element.")
                  warning(memo)
                  coverage <- coverage[1]
              }
              if(!is.numeric(coverage) && !is.null(coverage)){
                  memo <- paste("coverage is not numeric, attempting to coerce.")
                  warning(memo)
                  coverage <- as.numeric(coverage)
              }
              
              # status message
              if(verbose){
                  memo <- paste("Calculating a complex mutation burden.")
                  message(memo)
              }
              
              # obtain a data table of mutation counts on the sample level
              complexMutationCounts <- as.data.frame(table(primaryData[,c('sample', 'mutation')]))
              data.table::setDT(complexMutationCounts)
              
              # if coverage is not specified return just frequencies
              if(!is.numeric(coverage)){
                  if(verbose){
                      memo <- paste("coverage not specified, could not",
                                    "calculate the mutation burden... skipping")
                      message(memo)
                      complexMutationCounts$mutationBurden <- NA
                      return(complexMutationCounts)
                  }
              } 
              
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
#' @param genes character vector giving genes to keep
#' @param verbose Boolean for status updates
#' @return data.table object subset on gene if gene is not NULL. Entries are 
#' kept and or added if they are in the genes parameter.
#' @noRd
setMethod(f="geneSubset",
          signature="Waterfall",
          definition=function(object, genes, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData
              
              # Dont do anything if genes is null
              if(is.null(genes)) return(primaryData)
              
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

#' @rdname Waterfall-methods
#' @aliases mutHierarchySubset,Waterfall
#' @param object Object of class waterfall
#' @param mutationHierarchy character vector giving the order of mutations for
#' the hiearchy in order of most to least important
#' @param verbose Boolean for status updates
#' @return data.table object subset based on a mutation hierarchy keeping the
#' most important if there is more than one record for the same gene/sample.
#' @noRd
setMethod(f="mutHierarchySubset",
          signature="Waterfall",
          definition=function(object, verbose, ...){
              # grab the data
              primaryData <- object@primaryData
              mutationHierarchy <- object@MutationHierarchy
              
              # refactor the data frame
              primaryData$mutation<- factor(primaryData$mutation, levels=mutationHierarchy$mutation)
              
              # sort the data frame so that the duplicated call will remove the
              # proper mutation
              primaryData <- primaryData[order(primaryData$sample, primaryData$gene, primaryData$mutation),]
              
              # collapse the data on sample/gene
              primaryData <- primaryData[!duplicated(primaryData[, c("sample", "gene")]), ]
              
              # print status message
              if(verbose){
                  memo <- paste("Removed", nrow(object@primaryData)-nrow(primaryData),
                                "rows when setting the mutation hierarchy.")
                  message(memo)
              }
              
              return(primaryData)
          })

#' @rdname Waterfall-methods
#' @aliases recurrenceSubset,Waterfall
#' @param object Object of class waterfall
#' @param recurrence Numeric value specifying a recurrence cutoff to require,
#' genes not meeting this threshold are removed.
#' @param verbose Boolean for status updates
#' @return data.table object subset based on the recurrence of gene mutations.
#' @noRd
#' @importFrom stats na.omit
setMethod(f="recurrenceSubset",
signature="Waterfall",
definition=function(object, recurrence, verbose, ...){
    # access the part of the object we want to manipulate
    primaryData <- object@primaryData
    
    # Dont do anything if recurrence is null
    if(is.null(recurrence)) return(primaryData)
    
    # Perform quality checks
    if(!is.numeric(recurrence)){
        memo <- paste("argument supplied to recurrence is not of class",
                      "numeric, attempting to coerce!")
        warning(memo)
        recurrence <- as.numeric(as.character(recurrence))
    }
    if(length(recurrence) > 1){
        memo <- paste("argument supplied to recurrence has length > 1",
                      "only the first element was used.")
        warning(memo)
        recurrence <- recurrence[1]
    }
    
    # determine the frequency of gene mutations
    mutRecur <- primaryData[, count := .N, by = list(gene)]
    mutRecur <- unique(mutRecur[,c("gene", "count")])
    mutRecur <- stats::na.omit(mutRecur)
    mutRecur$prop <- mutRecur$count/nlevels(primaryData$sample)
    
    # If recurrence cutoff specified exceeds upper limit such that no
    # useful plot would be generated, reset recurrence cutoff
    maxRecur <- max(mutRecur$prop)
    if(maxRecur < recurrence){
        memo <- paste0("The recurrence cutoff specified exceeds the recurrence",
                       " seen in the data, resetting this value to equal max ",
                       "recurrence:", maxRecur)
        warning(memo)
        recurrence <- maxRecur
    }
    
    gene_above_recur <- mutRecur[mutRecur$prop >= recurrence,]$gene
    gene_below_recur <- mutRecur[mutRecur$prop < recurrence,]$gene
    
    # add NA to the end of 'gene_above_recurrence' vector, allowing for all
    # samples having NA as a gene name to be retained in the subset below
    gene_above_recur <- c(as.character(gene_above_recur), NA)
    
    # subset the original data frame based on the following: keep gene if it is
    # in the gene vector in "mutation_recurrence_subset"
    primaryData <- primaryData[(primaryData$gene %in% gene_above_recur), ]
    
    # print status message
    if(verbose){
        memo <- paste("Removing", length(unique(gene_below_recur)),
                      "genes not meeting the recurrence cutoff threshold.")
        message(memo)
    }
    
    return(primaryData)
})

#' @rdname Waterfall-methods
#' @aliases orderGenes,Waterfall
#' @param object Object of class waterfall
#' @param geneOrder Character vector specifying the order in which to plot
#' genes.
#' @param verbose Boolean for status updates
#' @return data.table object genes reordered
#' @noRd
setMethod(f="orderGenes",
          signature="Waterfall",
          definition=function(object, geneOrder, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData
             
              # print status message
              if(verbose){
                  memo <- paste("Setting gene order")
              }
              
              # order the genes based of frequency
              if(!is.null(geneOrder)) {
                  # perform quality checks on geneOrder
                  if(!is.character(geneOrder)){
                      memo <- paste("Argument supplied to gene order is not of",
                                    "class character, attempting to coerce.")
                      geneOrder <- as.character(geneOrder)
                      warning(memo)
                  }
                  if(any(duplicated(geneOrder))){
                      memo <- paste("Detected duplicated element in geneOrder,",
                                    "removing duplicates.")
                      geneOrder <- unique(geneOrder)
                      warning(memo)
                  }
                  
                  # if there are any genes in geneOrder not in x, remove those
                  gene_to_rmv <- geneOrder[!geneOrder %in% unique(primaryData$gene)]
                  if(length(gene_to_rmv) > 0){
                      memo <- paste("The following arguments to geneOrder were",
                                    "not found in the data:", toString(gene_to_rmv))
                      warning(memo)
                      geneOrder <- geneOrder[geneOrder %in% unique(primaryData$gene)]
                      primaryData$gene <- factor(primaryData$gene, levels=rev(geneOrder))
                      return(primaryData)
                  } else if(length(gene_to_rmv) == length(geneOrder)) {
                      memo <- paste0("Found no genes supplied to geneOrder in the",
                                     "data, check case.")
                      warning(memo)
                  }
              }
              
              # order based on mutation frequency
              gene_mutation_table <- table(primaryData[,c('gene', 'mutation')])
              geneOrder <- names(sort(rowSums(gene_mutation_table)))
              primaryData$gene <- factor(primaryData$gene, levels=geneOrder)
              return(primaryData)
          })

#' @rdname Waterfall-methods
#' @aliases maxGeneSubset,Waterfall
#' @param object Object of class waterfall
#' @param geneMax Integer specifying the maximum number of genes to be plotted.
#' @param verbose Boolean for status updates
#' @return data.table object subset to contain only the max number of genes.
#' @noRd
setMethod(f="maxGeneSubset",
          signature="Waterfall",
          definition=function(object, geneMax, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData

              # do nothing if null
              if(is.null(geneMax)){
                  return(primaryData)
              }
              
              # perform quality checks
              if(!is.numeric(geneMax)){
                  memo <- paste("geneMax is not numeric, attempting to convert",
                                "to an integer.")
                  geneMax <- as.integer(geneMax)
              }
              
              if(geneMax %% 1 != 0){
                  memo <- paste("geneMax is not a whole number, rounding.")
                  geneMax <- round(geneMax)
              }
              
              # limit to max genes
              keepGenes <- utils::tail(levels(primaryData$gene), geneMax)
              removeGenes <- unique(primaryData[!primaryData$gene %in% keepGenes,"gene"])
              primaryData <- primaryData[primaryData$gene %in% keepGenes,]
              
              # print status message
              if(verbose){
                  memo <- paste("geneMax is set to", geneMax, "removing",
                                nrow(removeGenes),"genes")
              }
              
              return(primaryData)
          })

#' @rdname Waterfall-methods
#' @aliases orderSamples,Waterfall
#' @param object Object of class waterfall
#' @param sampOrder Character vector specifying the order of samples
#' @param verbose Boolean for status updates
#' @return data.table object with samples reordered
#' @noRd
setMethod(f="orderSamples",
          signature="Waterfall",
          definition=function(object, sampleOrder, verbose, ...){
              # access the part of the object we want to manipulate
              primaryData <- object@primaryData
              
              # print status message
              if(verbose){
                  memo <- paste("setting the order of samples.")
                  message(memo)
              }
              
              # if sampleOrder is specified reorder on that
              if(!is.null(sampleOrder)){
                  # perform quality checks
                  if(!is.character(sampleOrder)){
                      memo <- paste("sampleOrder is not a character vector,",
                                    "attempting to coerce.")
                      warning(memo)
                      if(is.list(sampleOrder)){
                          sampleOrder <- unlist(sampleOrder)
                      }
                      sampleOrder <- as.character(sampleOrder)
                  }
                  if(sum(duplicated(sampleOrder)) != 0){
                      memo <- paste("Found duplicate elements in sampleOrder, uniquing!")
                      warning(memo)
                      sampleOrder <- unique(sampleOrder)
                  }
                  
                  # check if there are any samples specified not in the data
                  newSamples <- sampleOrder[!sampleOrder %in% levels(primaryData$sample)]
                  if(length(newSamples) != 0){
                      memo <- paste("The following samples were not detected",
                                    "in the data or its subsets:",
                                    toString(newSamples),
                                    ". Binding these to the data.")
                      warning(memo)
                      newSampleDT <- data.table::data.table("sample"=newSamples, "gene"=NA, "mutation"=NA, "label"=NA)
                      primaryData <- rbind(primaryData, newSampleDT, fill=TRUE)
                  }
                  
                  # status message
                  removeSample <- unique(primaryData$sample[!primaryData$sample %in% sampleOrder])
                  if(length(removeSample) != 0){
                      memo <- paste("Removing the samples:", toString(removeSample),
                                    "which were found in the data but not sampleOrder")
                      warning(memo)
                      primaryData <- primaryData[!primaryData$sample %in% removeSample,]
                  }
                  
                  # return the data with reordered samples
                  primaryData$sample <- factor(primaryData$sample, levels=sampleOrder)
                  return(primaryData)
              }
              
              # perform a hierarchical sort if sampleOrder is null
              # recast the data going from long format to wide format,
              #values in this data are counts of a mutation call
              wide_data <- reshape2::dcast(primaryData, sample ~ gene,
                                           fun.aggregate = length, value.var="mutation")
              
              # apply a boolean function to convert the data frame values to 1's and 0's
              values <- wide_data[,-1, drop=FALSE]
              sample <- wide_data[,1]
              values <- data.frame(apply(values, 2,
                                         function(x) as.numeric(as.logical(x))))
              wideBoolean <- cbind(sample, values)
              
              # reverse the columns so that genes with highest mutation's are 
              # listed first (assumes gene_sort has been run on the data frame)
              wideBoolean <- wideBoolean[,c(1, rev(2:ncol(wideBoolean)))]
              
              # if there are any NA values present in a sample at the gene put that
              # remove that sample and save (at this stage it should be samples)
              if(any(grepl("^NA.$", colnames(wideBoolean)))) {
                  # Find which column has the NA header
                  NA_index <- which(grepl("^NA.$", colnames(wideBoolean)))
                  
                  # Append NA column to end of data frame
                  NA_gene <- wide_boolean[,NA_index]
                  wide_boolean <- wide_boolean[,-NA_index]
                  
                  # Save copy and remove samples with no mutations,
                  # these will be added to the end
                  samp_no_mut <- wideBoolean[rowSums(wideBoolean[2:ncol(wideBoolean)]) == 0,]$sample
                  samp_no_mut <- as.character(samp_no_mut)
                  wideBoolean <- wideBoolean[!wideBoolean$sample %in% samp_no_mut,]
              } else {
                  samp_no_mut <- NULL
              }
              
              # hiearchial sort on all column's (i.e. genes) such that samples are
              # rearranged if there is a mutation in that gene
              sampleOrder <- wideBoolean[do.call(order, as.list(-wideBoolean[2:ncol(wideBoolean)])),]$sample
              
              # Put those samples not in sample order in from the original levels of the
              # data (these are samples with no mutations)
              not_in <- as.character(levels(sampleOrder)[which(!levels(sampleOrder) %in% sampleOrder)])
              not_in <- not_in[!not_in %in% samp_no_mut]
              sampleOrder <- c(as.character(sampleOrder), as.character(not_in))
              
              # Put those samples with no mutations back in
              if(!is.null(samp_no_mut)){
                  sampleOrder <- c(sampleOrder, samp_no_mut)
              }
              
              # set the new sample order
              primaryData$sample <- factor(primaryData$sample, levels=sampleOrder)
              return(primaryData)
          })

#' @rdname Waterfall-methods
#' @aliases buildMutationPlot,Waterfall
#' @param object Object of class waterfall
#' @param plotA String specifying the type of plot for the top sub-plot, one of
#' "burden", "frequency", or NULL for a mutation burden (requires coverage to be
#' specified), frequency of mutations, or no plot respectively.
#' @param plotATally String specifying one of "simple" or "complex" for a
#' simplified or complex tally of mutations respectively.
#' @param plotALayers list of ggplot2 layers to be passed to the plot.
#' @param verbose Boolean for status updates
#' @return gtable object containing the top sub-plot.
#' @noRd
#' @import ggplot2
#' @importFrom gtable gtable
setMethod(f="buildMutationPlot",
          signature="Waterfall",
          definition=function(object, plotA, plotATally, plotALayers, verbose, ...){
              # grab only the first element for parameters
              plotA <- plotA[1]
              plotATally <- plotATally[1]
              
              # if plot type is null return an empty gtable
              if(is.null(plotA)) return(gtable::gtable())

              # print status message
              if(verbose) {
                  memo <- paste("Constructing top sub-plot")
                  message(memo)
              }
              # extract the data for the type of plot we need
              if(toupper(plotATally) == toupper("simple")){
                  mutationData <- object@simpleMutationCounts
              } else if(toupper(plotATally) == toupper("complex")) {
                  mutationData <- object@complexMutationCounts
              }
              
              # perform quality checks
              if(!is.null(plotALayers) && !is.list(plotALayers)){
                  memo <- paste("plotALayers is not a list... attempting to coerce.")
                  warning(memo)
                  plotALayers <- as.list(plotALayers)
                  if(any(lapply(plotALayers, function(x) ggplot2::is.ggproto(x) | ggplot2::is.theme(x)))){
                      memo <- paste("plotALayers is not a list of ggproto or ",
                                    "theme objects... setting plotALayers to NULL")
                      warning(memo)
                      plotALayers <- NULL
                  }
              }
              if(!is.null(plotA) && !toupper(plotA) %in% toupper(c("frequency", "burden"))){
                  memo <- paste("plotA is not set to \"frequency\", \"burden\"",
                                " or NULL... defaulting to \"frequency\".")
                  message(memo)
                  plotA <- "frequency"
              }
              if(all(is.na(mutationData$mutationBurden)) && toupper(plotA) == toupper("burden")){
                  memo <- paste("plotA is set to:",toString(plotA),"but could",
                                "not find calculated mutation burden, please",
                                "specify coverage!, Resetting plotA to frequency.")
                  warning(memo)
                  plotA <- "frequency"
              }
              if(!toupper(plotATally) %in% toupper(c("simple", "complex"))){
                  memo <- paste("plotATally is not set to either \"simple\" or",
                                " \"complex\"... defaulting to \"simple\"")
                  message(memo)
                  plotATally <- "simple"
              }
              
              # make sure sample levels match primaryData for plotting
              mutationData <- mutationData[mutationData$sample %in% unique(object@primaryData$sample),]
              mutationData$sample <- factor(mutationData$sample, levels=levels(object@primaryData$sample))
              if(toupper(plotATally) == toupper("complex")) {
                  mutationData$mutation <- factor(mutationData$mutation, levels=levels(object@primaryData$mutation))
              } else if(toupper(plotATally) == toupper("simple")) {
                  mutationData$mutation <- factor(mutationData$mutation, levels=c("Non Synonymous","Synonymous"))
              }
              
              ############# set ggplot2 layers #################################
              
              # theme
              plotTheme <- theme(axis.ticks.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.title.x = element_blank(),
                                 legend.title = element_text(size=14),
                                 axis.text.y = element_text(colour = "black"),
                                 axis.title.y = element_text(colour = "black"),
                                 panel.background = element_blank(),
                                 # panel.grid.minor.y = element_line(colour = "black"),
                                 panel.grid.major.y = element_line(colour = "grey80"),
                                 panel.grid.minor.x = element_blank(),
                                 panel.grid.major.x = element_blank(),
                                 panel.border = element_rect(fill = NA)
              )
              # legend
              if(toupper(plotATally) == toupper("simple")){
                  plotLegend <- scale_fill_manual(name="Translational Effect",
                                                  values=c("Synonymous"="red", "Non Synonymous"="blue"),
                                                  breaks=c("Synonymous", "Non Synonymous"),
                                                  drop=FALSE)
              } else if(toupper(plotATally) == toupper("complex")){
                  plotLegend <- scale_fill_manual(name="Translational Effect",
                                                  values=object@MutationHierarchy$color,
                                                  breaks=object@MutationHierarchy$mutation,
                                                  drop=FALSE)
              }
              
              # titles
              if(toupper(plotA) == toupper("frequency")) {
                  plotTitleY <- ylab("Mutation Frequency")
              } else if(toupper(plotA) == toupper("burden")) {
                  plotTitleY <- ylab("Mutations\nper MB")
              }
              
              #  geom definition
              plotGeom <- geom_bar(stat='identity', alpha=.75, width=1)
              
              # plot
              if(toupper(plotA) == toupper("frequency")) {
                  mutPlot <- ggplot(mutationData, aes_string(x='sample', y='Freq', fill='mutation'))
              } else if(toupper(plotA) == toupper("burden")) {
                  mutPlot <- ggplot(mutationData, aes_string(x='sample', y='mutationBurden', fill='mutation'))
              }
              mutPlot <- mutPlot + plotGeom + plotTitleY + plotLegend + plotTheme + plotALayers
              
              # convert to gtable grob
              plotGrob <- ggplotGrob(mutPlot)
              return(plotGrob)
          })

#' @rdname Waterfall-methods
#' @aliases constructGeneData,Waterfall
#' @param object Object of class waterfall
#' @param verbose Boolean for status updates
#' @return data.table object containing summarized gene level data
#' @noRd
setMethod(f="constructGeneData",
          signature="Waterfall",
          definition=function(object, verbose, ...){
              # extract the data to work with
              geneData <- object@primaryData
              
              # status message
              if(verbose){
                  memo <- paste("Constructing GeneData")
              }
              
              # construct geneData
              geneData <- geneData[, count := .N, by = list(gene, mutation)]
              geneData <- unique(geneData[,c("gene", "mutation", "count")])
              geneData$gene <- factor(geneData$gene, levels=levels(object@primaryData$gene))
              
              return(geneData)
          })

#' @rdname Waterfall-methods
#' @aliases buildGenePlot,Waterfall
#' @param object Object of class waterfall
#' @param plotB String specifying the type of plot for the left sub-plot, one of
#' "proportion", "frequency", or NULL for a plot of gene proportions frequencies
#' , or no plot respectively.
#' @param plotBTally String specifying one of "simple" or "complex" for a
#' simplified or complex tally of genes respectively.
#' @param plotBLayers list of ggplot2 layers to be passed to the plot.
#' @param verbose Boolean for status updates
#' @return gtable object containing the left sub-plot.
#' @noRd
#' @import ggplot2
#' @importFrom gtable gtable
setMethod(f="buildGenePlot",
          signature="Waterfall",
          definition=function(object, plotB, plotBTally, plotBLayers, verbose, ...){
              # grab only the first element for parameters
              plotB <- plotB[1]
              plotBTally <- plotBTally[1]
              
              # if plot type is null return an empty gtable
              if(is.null(plotB)) return(gtable::gtable())
              
              # extract the data needed for this plot
              geneData <- object@geneData
              
              # print status message
              if(verbose) {
                  memo <- paste("Constructing left sub-plot")
                  message(memo)
              }
              browser()
              # perform quality checks
              if(!is.null(plotBLayers) && !is.list(plotBLayers)){
                  memo <- paste("plotBLayers is not a list... attempting to coerce.")
                  warning(memo)
                  plotBLayers <- as.list(plotBLayers)
                  if(any(lapply(plotBLayers, function(x) ggplot2::is.ggproto(x) | ggplot2::is.theme(x)))){
                      memo <- paste("plotBLayers is not a list of ggproto or ",
                                    "theme objects... setting plotBLayers to NULL")
                      warning(memo)
                      plotBLayers <- NULL
                  }
              }
              if(!is.null(plotB) && !toupper(plotB) %in% toupper(c("frequency", "proportion"))){
                  memo <- paste("plotB is not set to \"frequency\", \"proportion\"",
                                " or NULL... defaulting to \"proportion\".")
                  message(memo)
                  plotB <- "proportion"
              }
              if(!toupper(plotBTally) %in% toupper(c("simple", "complex"))){
                  memo <- paste("plotBTally is not set to either \"simple\" or",
                                " \"complex\"... defaulting to \"simple\"")
                  message(memo)
                  plotBTally <- "simple"
              }
              
              # determine proportion of gene mutations
              sampleCount <- nlevels(object@primaryData$sample)
              geneData$proportion <- geneData$count/sampleCount * 100
              
              ################ ggplot2 #########################################
              
              # theme
              plotTheme <- theme(axis.text.y=element_text(colour='black', face='italic'),
                                 axis.title.y=element_blank(),
                                 legend.position=('none'))
              # titles
              if(plotB == "frequency"){
                  plotYlabel <- ylab('# Mutant')
              }else if(plotB == "proportion"){
                  plotYlabel <- ylab('% Mutant')
              }
              
              # legend
              if(plotBTally == "simple"){
                  plotLegend <- scale_fill_manual(name="Translational Effect",
                                                  values=object@MutationHierarchy$color,
                                                  breaks=object@MutationHierarchy$mutation,
                                                  drop=FALSE)
              }else if(plotBTally == "complex"){
                  plotLegend <- scale_fill_manual(name="Translational Effect",
                                                  values=object@MutationHierarchy$color,
                                                  breaks=object@MutationHierarchy$mutation,
                                                  drop=FALSE)
              }
              
              # geom definition
              plotGeom <- geom_bar(position='stack', alpha=.75, width=1, stat='identity')
              
              # plot
              if(plotB == "proportion"){
                  genePlot <- ggplot(geneData, aes_string(x='gene', y='proportion', fill='mutation'))
              }else if(plotB == "frequency"){
                  genePlot <- ggplot(geneData, aes_string(x='gene', y='count', fill='mutation'))
              }
              genePlot <- genePlot + plotGeom + theme_bw() + coord_flip() + plotTheme +
                  plotYlabel + scale_y_reverse() + plotLegend + plotBLayers
              
              plotGrob <- ggplotGrob(genePlot)
          })

#' @rdname formatClinicalData-methods
#' @aliases formatClinicalData,MutationAnnotationFormat
#' @noRd
#' @importFrom data.table setDT
#' @importFrom data.table is.data.table
#' @importFrom data.table melt
setMethod(f="formatClinicalData",
          signature="Waterfall",
          definition=function(object, verbose, ...){
              
              # extract the data we need
              clinicalData <- object@ClinicalData
              primaryData <- object@primaryData
              
              # print status message
              if(verbose){
                  memo <- paste("Formatting clinical data")
                  message(memo)
              }
              
              # remove clinical samples not found in primary data
              clinicalData <- clinicalData[clinicalData$sample %in% primaryData$sample,]
              
              # check no samples were removed
              removedSamp <- unique(object@ClinicalData$sample[!object@ClinicalData$sample %in% clinicalData$sample])
              memo <- paste("Removed", length(removedSamp), "samples from the clinical data",
                            "not found in the primary waterfall data! These samples are:",
                            toString(removedSamp))
              if(removedSamp != 0){
                  warning(memo)
              } else if(verbose){
                  message(memo)
              }
              
              # set levels of clinicalData to match primaryData
              clinicalData$sample <- factor(clinicalData$sample, levels=levels(primaryData$sample))
              
              # return the formated data
              return(clinicalData)
          })

#' @rdname Waterfall-methods
#' @aliases buildWaterfallPlot,Waterfall
#' @param object Object of class waterfall
#' @param verbose Boolean for status updates
#' @param gridOverlay Boolean specifying if a grid should be overlayed on the
#' waterfall plot.
#' @param drop Boolean specifying if unused mutations should be dropped from the
#' legend.
#' @return gtable object containing the main plot.
#' @noRd
#' @import ggplot2
#' @importFrom gtable gtable
setMethod(f="buildWaterfallPlot",
          signature="Waterfall",
          definition=function(object, gridOverlay, drop, labelSize, labelAngle, sampleNames, xTitle, verbose, ...){
              # extract the data we need
              primaryData <- object@primaryData
              paletteData <- object@MutationHierarchy
              
              # print status message
              if(verbose){
                  memo <- paste("Building the main plot")
                  message(memo)
              }
              
              ######### start building the plot ################################
              # grid overlay
              if(gridOverlay){
                  verticalPlotGrid <- geom_vline(xintercept = seq(.5, nlevels(primaryData$sample),
                                                                  by=1),
                                                 linetype='solid', colour='grey80', size=.01)
                  
                  if(length(unique(primaryData$gene)) == 1)
                  {
                      horizontalPlotGrid <- geom_blank()
                  } else {
                      horizontalPlotGrid <- geom_hline(yintercept = seq(1.5, length(unique(primaryData$gene)),
                                                                        by=1),
                                                       linetype='solid', colour='grey80',
                                                       size=.01)
                  }
                  plotGridOverlay <- list(horizontalPlotGrid, verticalPlotGrid)
              } else {
                  plotGridOverlay <- geom_blank()
              }
              
              # construct legend
              if(drop)
              {
                  valueMap <- paletteData$color
                  names(valueMap) <- paletteData$mutation
                  plotLegend <- scale_fill_manual(name="Mutation Type",
                                                  values=valueMap,
                                                  breaks=paletteData$mutation,
                                                  labels=paletteData$label,
                                                  drop=TRUE)
              } else {
                  valueMap <- paletteData$color
                  names(valueMap) <- paletteData$mutation
                  plotLegend <- scale_fill_manual(name="Mutation Type",
                                                  values=valueMap,
                                                  breaks=paletteData$mutation,
                                                  labels=paletteData$label,
                                                  drop=FALSE)
              }
              
              # plot titles
              plotXLabel <- xlab(paste0('Sample (n=', nlevels(primaryData$sample), ')'))
              
              if(all(is.na(primaryData$label)))
              {
                  label <- geom_text(data=primaryData,
                                     mapping=aes_string(x='sample', y='gene',
                                                        label='label'),
                                     size=labelSize, colour='white',
                                     angle=labelAngle)
              } else {
                  label <- geom_blank()
              }
              
              # base theme
              plotTheme <- theme(axis.ticks=element_blank(),
                                 panel.grid.major=element_blank(),
                                 panel.grid.minor=element_blank(),
                                 axis.title.y=element_blank(),
                                 panel.background=element_rect(fill='white',
                                                               colour='white'),
                                 axis.text.y=element_blank(),
                                 plot.title=element_blank(),
                                 legend.title=element_text(size=14))
              
              if(sampleNames){
                  sampleLabels <- theme(axis.text.x=element_text(angle=50, hjust=1))
              } else {
                  sampleLabels <- theme(axis.text.x=element_blank(),
                                        axis.title.x=element_blank())
              }
              if(xTitle){
                  plotXtitle <- theme(axis.title.x=element_text(size=20))
              } else {
                  plotXtitle <- theme(axis.title.x=element_blank())
              }
              
              # define geom
              plotGeom <- geom_tile(aes_string(fill='mutation'), position="identity")
              
              # define plot
              waterfallPlot <- ggplot(primaryData, aes_string('sample', 'gene'))
              
              # combine plot elements
              waterfallPlot <- waterfallPlot + plotGeom + plotGridOverlay +
                  plotLegend + plotXLabel + plotTheme + sampleLabels + plotXtitle
              
              # covert to grob
              waterfallGrob <- ggplotGrob(waterfallPlot)
              return(waterfallGrob)
          })

#' @rdname buildClinicalPlot-methods
#' @aliases buildClinicalPlot,Waterfall
#' @param clinicalLayers List of ggplot2 layers to add to the clinical plot.
#' @noRd
#' @import ggplot2
#' @importFrom gtable gtable
setMethod(f="buildClinicalPlot",
          signature="Waterfall",
          definition=function(object, clinicalLayers, ...){
              # extract necessary data
              clinicalData <- object@ClinicalData
              
              # if clinical data is empty return empty gtable
              if(nrow(clinicalData) == 0){
                  return(gtable::gtable())
              }
              
              # construct a plot
              clinicalXLabel <- xlab(paste0("Sample n=", length(unique(clinicalData$sample))))
              clinicalPlot <- ggplot(clinicalData, aes_string(x='sample',
                                                              y='variable',
                                                              fill='value')) +
                  clinicalLayers + clinicalXLabel + theme(axis.text.x=element_blank())
              
              # contruct grob
              clinicalGrob <- ggplotGrob(clinicalPlot)
          })

#' @rdname arrangeWaterfallPlot-methods
#' @aliases arrangeWaterfallPlot,Waterfall
#' @param sectionHeights Relative heights of each plot section (should sum to one).
#' @param sectionWidths Relative widths of each plot section (should sum to one).
#' @noRd
#' @importFrom grid nullGrob
#' @importFrom gridExtra arrangeGrob
setMethod(f="arrangeWaterfallPlot",
          signature="Waterfall",
          definition=function(object, sectionHeights, sectionWidths, verbose, ...){
              
              # grab the data we need
              plotA <- object@PlotA
              plotB <- object@PlotB
              plotC <- object@PlotC
              plotD <- object@PlotD
              
              # set default widths
              if(is.null(sectionWidths)){
                  sectionWidths <- c(.25, .75)
              } else if(length(sectionWidths) != 2){
                  memo <- paste("sectionWidths should be of length 2... Using",
                                "default values!")
                  warning(memo)
                  sectionWidths <- c(.25, .75)
              } else if(!all(is.numeric(sectionWidths))) {
                  memo <- paste("sectionWidths must be numeric... Using",
                                "default values!")
                  warning(memo)
                  sectionWidths <- c(.25, .75)
              }
              
              # set up a blank plot for alignment purposes
              emptyPlot <- grid::nullGrob()
              
              ## Strip out legends and plot separately
              ## https://github.com/baptiste/gridextra/wiki/arranging-ggplot#legends
              #ind_legend <- grep("guide", plotC$layout$name)
              #plotC_legend <- plotC[["grobs"]][[ind_legend]]
              #plotC_width <- sum(plotC_legend$width)
              #heatmap_grob <- ggplot2::ggplotGrob(grid.draw(plotC) + theme(legend.position="none"))
              
              # align the heights of plotB and plotC if they exist
              if(length(plotB) != 0){
                  maxHeight <- grid::unit.pmax(plotB$heights, plotC$heights)
                  plotB$heights <- as.list(maxHeight)
                  plotC$heights <- as.list(maxHeight)
              }
              
              # Obtain the max width for relevant plots
              plotList4Width <- list(plotA, plotC, plotD)
              plotList4Width <- plotList4Width[lapply(plotList4Width, length) > 0]
              plotList4Width <- lapply(plotList4Width, function(x) x$widths)
              maxWidth <- do.call(grid::unit.pmax, plotList4Width)
              
              # section plots into rows rows, also define preliminary heights
              plotList <- list()
              plotHeight <- list()
              # set width for first row of plots
              if(length(plotA) != 0){
                  plotA$widths <- as.list(maxWidth)
                  firstRowPlots <- gridExtra::arrangeGrob(emptyPlot, plotA, ncol=2, widths=sectionWidths)
                  plotList[["firstRow"]] <- firstRowPlots
                  plotHeight[["firstRow"]] <- .2
              }
              
              # set widths for second row of plots
              plotC$widths <- as.list(maxWidth)
              if(length(plotB) != 0){
                  secondRowPlots <- gridExtra::arrangeGrob(plotB, plotC, ncol=2, widths=sectionWidths)
              } else {
                  secondRowPlots <- gridExtra::arrangeGrob(emptyPlot, plotC, ncol=2, widths=sectionWidths)
              }
              plotList[["secondRow"]] <- secondRowPlots
              plotHeight[["secondRow"]] <- .8
              
              # set widths for third row of plots
              if(length(plotD) != 0){
                  plotD$widths <- as.list(maxWidth)
                  lastRowPlots <- gridExtra::arrangeGrob(emptyPlot, plotD, ncol=2, widths=sectionWidths)
                  plotList[["lastRow"]] <- lastRowPlots
                  plotHeight[["lastRow"]] <- .2
              }
              
              # set section heights based upon the number of sections
              plotHeight <- as.vector(do.call(rbind, plotHeight))
              plotHeight[which(plotHeight == min(plotHeight))] <- min(plotHeight)/length(which(plotHeight == min(plotHeight)))
              
              if(is.null(sectionHeights)){
                  sectionHeights <- plotHeight
              } else if(length(sectionHeights) != length(plotList)){
                  memo <- paste("There are", length(sectionHeights), "provided",
                                "but", length(plotList), "vertical sections...",
                                "using default values!")
                  warning(memo)
                  sectionHeights <- plotHeight
              } else if(!all(is.numeric(sectionHeights))) {
                  memo <- paste("sectionHeights must be numeric... Using",
                                "default values!")
                  warning(memo)
                  sectionHeights <- plotHeight
              }
              
              # arrange the final plot
              # http://stackoverflow.com/questions/6681145/how-can-i-arrange-an-arbitrary-number-of-ggplots-using-grid-arrange
              finalPlot <- do.call(gridExtra::arrangeGrob, c(plotList, list(ncol=1, heights=sectionHeights)))
              
              return(finalPlot)
          })

#' @rdname drawPlot-methods
#' @aliases drawPlot,Waterfall
#' @importFrom grid grid.draw
setMethod(
    f="drawPlot",
    signature="Waterfall",
    definition=function(object, ...){
        mainPlot <- object@Grob
        grid::grid.draw(mainPlot)
    }
)