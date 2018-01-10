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
#' @slot primaryData data.table object holding the primary data plotted
#' @slot geneData data.table object holding transcript and domain coordinates
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
LolliplotData <- function(object, gene, transcript, species, host, verbose){
    
    # convert object to Lolliplot format
    lolliplotData <- toLolliplot(object, verbose=verbose)
    
    # filter the data by gene as a first step if no transcript is imediately available
    # this will reduce runtime by dramatically subsetting the data almost immediately
    lolliplotData <- filterByGene(lolliplotData, transcript=transcript, gene=gene, verbose=verbose)
    
    # annotate data with ensembl gene if necessary
    lolliplotData <- annotateGene(lolliplotData, species=species, host=host, verbose=verbose)
    
    # annotate data with ensembl transcripts if necessary
    lolliplotData <- annotateTranscript(lolliplotData, species=species, host=host, verbose=verbose)
    browser()
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
              index <- object$gene %in% gene
              object <- object[index,]
              object$gene <- factor(object$gene, levels=unique(object$gene))
              
              return(object)
          })

#' @rdname annotateTranscript-methods
#' @aliases annotateTranscript
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @importFrom biomaRt useMart
#' @importFrom biomaRt listDatasets
#' @importFrom biomaRt useDataset
#' @importFrom biomaRt getBM
#' @noRd
setMethod(f="annotateTranscript",
          signature="data.table",
          definition=function(object, species, host, verbose, ...){
              
              # first check if a transcript column exists
              if("transcript" %in% colnames(object)){
                  
                  # next check to see if they appear to be ensembl transcripts
                  if(any(grepl("ENST00", object$transcript))){
                      
                      if(verbose){
                          memo <- paste("Found a transcript column with ensembl transcript annotations",
                                        "skipping transcript annotation")
                          message(memo)
                      }
                      object$ensembl_transcript_id <- object$transcript
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
              
              # select proper data set given regexp print warnings if unexpected out occur
              dataset <- biomaRt::listDatasets(ensembl_mart)$dataset
              index <- which(grepl(species, dataset))
              if(length(index)>1)
              {
                  memo <- paste(toString(species), " Matches more than one dataset for the",
                                 " ensembl mart, please specify a species in the, ",
                                 "following format: hsapiens")
                  stop(memo)
              } else if(length(index)==0) {
                  memo <- paste(toString(species), " does not appear to be supported by biomaRt",
                                 "if you beleive this to be in error please modify", 
                                 "you're input to to conform to this format: hsapiens")
                  stop(memo)
              }
              
              ensembl_mart <- biomaRt::useDataset(as.character(dataset[index]),
                                                  mart=ensembl_mart)
             
              # Apply various filters using vector of values
              # biomaRt::listFilters(ensembl_mart)
              filters <- c("ensembl_gene_id")
              values <- c(as.character(object$ensembl_gene_id))
              
              # Select attributes to retrieve (protein domain, start, stop)
              # biomaRt::listAttributes(ensembl_mart)
              attributes <- c("ensembl_transcript_id")
             
              # Retrieve data
              result <- biomaRt::getBM(attributes=attributes, filters=filters,
                                       values=values, mart=ensembl_mart)
              
              # add in the data and return the structure
              object$ensembl_transcript_id <- result
              
              return(object)
          })

#' @rdname annotateGene-methods
#' @aliases annotateGene
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @importFrom biomaRt useMart
#' @importFrom biomaRt listDatasets
#' @importFrom biomaRt useDataset
#' @importFrom biomaRt getBM
#' @noRd
setMethod(f="annotateGene",
          signature="data.table",
          definition=function(object, species, host, verbose, ...){
              
              # next check to see if they appear to be ensembl transcripts
              if(all(grepl("ENSG00", object$gene))){
                  
                  if(verbose){
                      memo <- paste("Found a gene column with ensembl gene annotations",
                                    "skipping gene annotation")
                      message(memo)
                  }
                  object$ensembl_gene_id <- object$gene
                  return(object)
              }

              # annotate the gene using the gene
              if(verbose){
                  memo <- paste("Annotating ensembl genes via biomaRt")
              }
              
              # load the appropriate mart
              ensembl_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                                               host=host)
              
              # select proper data set given regexp print warnings if unexpected out occur
              dataset <- biomaRt::listDatasets(ensembl_mart)$dataset
              index <- which(grepl(species, dataset))
              if(length(index)>1)
              {
                  memo <- paste(toString(species), " Matches more than one dataset for the",
                                " ensembl mart, please specify a species in the, ",
                                "following format: hsapiens")
                  stop(memo)
              } else if(length(index)==0) {
                  memo <- paste(toString(species), " does not appear to be supported by biomaRt",
                                "if you beleive this to be in error please modify", 
                                "you're input to to conform to this format: hsapiens")
                  stop(memo)
              }
              
              ensembl_mart <- biomaRt::useDataset(as.character(dataset[index]),
                                                  mart=ensembl_mart)
              
              # Apply various filters using vector of values
              # biomaRt::listFilters(ensembl_mart)
              filters <- c("hgnc_symbol")
              values <- c(as.character(object$gene))
              
              # Select attributes to retrieve (protein domain, start, stop)
              # biomaRt::listAttributes(ensembl_mart)
              attributes <- c("ensembl_gene_id")
              
              # Retrieve data
              result <- biomaRt::getBM(attributes=attributes, filters=filters,
                                       values=values, mart=ensembl_mart)
              
              # add in the data and return the structure
              object$ensembl_gene_id <- result
              
              return(object)
          })