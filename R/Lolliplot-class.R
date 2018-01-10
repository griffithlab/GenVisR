################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Class Lolliplot
#' 
#' An S4 class for the lolliplot object
#' @name Lolliplot-class
#' @rdname Lolliplot-class
#' @slot PlotA gtable object for the top sub-plot
#' @slot PlotB gtable object for the bottom sub-plot
#' @slot Grob gtable object storing the arranged plot
#' @slot primaryData data.table object storing the primary data
#' @slot geneData data.table object storing gene and domain coordinates
#' @exportClass Lolliplot
#' @import methods
#' @importFrom gtable gtable
#' @importFrom data.table data.table
methods::setOldClass("gtable")
setClass("Lolliplot",
         representation=representation(PlotA="gtable",
                                       PlotB="gtable",
                                       Grob="gtable",
                                       primaryData="data.table",
                                       geneData="data.table"),
         validity=function(object){
             
         })

#' Constructor for the Lolliplot class.
#' 
#' @name Lolliplot
#' @rdname Lolliplot-class
#' @param input Object of class MutationAnnotationFormat, GMS, VEP, or a data.table with appropriate columns
#' @param gene Character string specifying a gene to plot
#' @param transcript Character string specifying the ensembl transcript for which to plot, should be a transcript which corresponds
#' to the gene parameter.
#' @param species Character string specifying a species when using biomaRt queries
#' @param host Character string specifying a host to connect to when using biomaRt queries
#' @export
Lolliplot <- function(input, gene, transcript, species="hsapiens", host="www.ensembl.org", verbose=FALSE){
    
    # Obtain and format the data
    data <- LolliplotData(input, gene=gene, transcript=transcript, species=species, host=host, verbose=verbose)
    
    
}

#!!!!!!!!!!!!!!!!!!!!!!!! Private Classes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Private Class LolliplotData
#' 
#' An S4 class for the data of the Lolliplot object
#' @name LolliplotData-class
#' @rdname LolliplotData-class
#' @slot primaryData
#' @slot geneData
#' @import methods
#' @importFrom data.table data.table
#' @noRd
setClass("LolliplotData",
         representation=representation(primaryData="data.table",
                                       geneData="data.table"),
         validity=function(object){
             
         })

#' Constructor for the LolliplotData class
#' 
#' @name LolliplotData
#' @rdname LolliplotData-class
#' @param object Object of class MutationAnnotationFormat, GMS, VEP, or a data.table with appropriate fields
#' @noRd
LolliplotData <- function(object, gene, species, host, verbose){
    
    # convert object to Lolliplot format
    lolliplotData <- toLolliplot(object, verbose)
    
    # filter the data by gene as a first step if no transcript is imediately available
    # this will reduce runtime by dramatically subsetting the data almost immediately
    lolliplotData <- filterByGene(lolliplotData, transcript, gene, verbose)
    browser()
    # annotate data with ensembl transcripts if necessary
    lolliplotData <- annotateTranscript(lolliplotData, species, host, verbose)
    
    # filter data to only one transcript
    #lolliplotData <- filterByTranscript(lolliplotData, transcript, verbose)
    
    # grab the transcript coordinates
    
    # set the initial mutation heights (tier 1)
    
    # set the mutations for the second tier
}

#' PrivateClass LolliplotPlots
#' 
#' An S4 class for the grobs of the Lolliplot object
#' @name LolliplotPlots-class
#' @rdname LolliplotPlots-class
#' @slot PlotA gtable object for the top plot.
#' @slot PlotB gtable object for the bottom plot.
#' @import methods
#' @importFrom gtable gtable
#' @noRd
setClass("LolliplotPlots", 
         representation=representation(PlotA="gtable",
                                       PlotB="gtable"),
         validity=function(object){
             
         })

#' Constructor for the LolliplotPlots class
#' 
#' @name LolliplotPlots
#' @rdname LolliplotPlots-class
#' @param object Object of class LolliplotData
#' @noRd
LolliplotPlots <- function(object, verbose){
    
}

################################################################################
########################### Accessor function definitions ######################

################################################################################
########################### Method function definitions ########################

#' @rdname filterByGene-methods
#' @aliases filterByGene
#' @param object Object of class data.table
#' @param gene Character string specifying a gene name
#' @param verbose Boolean specifying if status messages should be reported
#' @noRd
setMethod(f="filterByGene",
          signature="data.table",
          definition=function(object, gene, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Filtering data by the gene:", toString(gene))
                  message(memo)
              }
              
              # perform quality checks to the parameter gene
              if(!is.character(gene)){
                  memo <- paste("argument supplied to the parameter gene is not",
                                "a character vector... attempting to coerce")
                  gene <- as.character(gene)
                  warning(memo)
              }
              if(length(gene) != 1){
                  memo <- paste("argument supplied to the parameter gene has a",
                                "length greater than 1... using only the first element")
                  gene <- gene[1]
                  warning(memo)
              }
              
              # check that the gene to filter is in the data
              if(!any(object$gene %in% gene)){
                  memo <- paste("The string given to the parameter gene does not",
                                "exist in the data!")
                  stop(memo)
              }
              
              # perform the filtering and refactor the column
              object <- object[object$gene %in% gene,]
              object$gene <- factor(object$gene, levels=unqiue(object$gene))
              
              return(object)
          })

#' @rdname annotateTranscript-methods
#' @aliases annotateTranscript
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @noRd
setMethod(f="annotateTranscript",
          signature="data.table",
          definition=function(object, species, host, verbose, ...){
              
              # first check if a transcript column exists
              if("transcript" %in% colnames(object)){
                  
                  # next check to see if they appear to be ensembl transcripts
                  if(any(grepl("ENST", object$transcript))){
                      
                      if(verbose){
                          memo <- paste("Found a transcript column with ensembl transcript annotations",
                                        "skipping transcript annotation")
                          message(memo)
                      }
                      return(object)
                  }
              }
              
              # annotate the transcripts using the gene
              if(verbose){
                  memo <- paste("Annotating ensembl transcripts via biomaRt")
              }
              
              # load the appropriate mart
              ensembl_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                                               host=host)
              browser()
              
              return(object)
          })