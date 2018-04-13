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
#' @param transcript Character string specifying the ensembl transcript for which to plot, should be a transcript which corresponds
#' to the gene parameter.
#' @param species Character string specifying a species when using biomaRt queries
#' @param host Character string specifying a host to connect to when using biomaRt queries
#' @param txdb A bioconoductor txdb object to annotate amino acid positions, required only if amino acid changes are missing (see details).
#' @param BSgenome A bioconductor BSgenome object to annotate amino acid positions, required only if amino acid changes are missing (see details).
#' @param emphasize Character vector specifying a list of mutations to emphasize.
#' @param DomainPalette Character vector specifying the colors used for encoding protein domains
#' @param MutationPalette Character vector specifying the colors used for encoding mutations
#' @param plotALayers list of ggplot2 layers to be passed to the density plot.
#' @param plotBLayers list of ggplot2 layers to be passed to the lolliplot.
#' @param sectionHeights Numeric vector specifying relative heights of each plot section,
#' should sum to one. Expects a value for each section.
#' @param labelAA Boolean specifying if labels should be added to emphasized mutations
#' @param verbose Boolean specifying if status messages should be reported.
#' @export
Lolliplot <- function(input, transcript=NULL, species="hsapiens", host="www.ensembl.org", txdb=NULL, BSgenome=NULL, emphasize=NULL, DomainPalette=NULL, MutationPalette=NULL, labelAA=TRUE, plotALayers=NULL, plotBLayers=NULL, sectionHeights=NULL, verbose=FALSE){
    
    # Obtain and format the data
    data <- LolliplotData(input, transcript=transcript, species=species, host=host, txdb=txdb, BSgenome=BSgenome, emphasize=emphasize, verbose=verbose)
    
    # construct the plots
    lolliplotPlots <- LolliplotPlots(data, labelAA=labelAA, DomainPalette=DomainPalette, MutationPalette=MutationPalette, plotALayers=plotALayers, plotBLayers, verbose=verbose)
    
    # align the plots
    lolliplotGrob <- arrangeLolliplotPlot(lolliplotPlots, sectionHeights=sectionHeights, verbose=verbose)
    
    new("Lolliplot", PlotA=getGrob(lolliplotPlots, index=1), PlotB=getGrob(lolliplotPlots, index=2), Grob=lolliplotGrob, primaryData=getData(data, name="primaryData"), geneData=getData(data, name="geneData"))
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
LolliplotData <- function(object, transcript, species, host, txdb, BSgenome, emphasize, verbose){
    
    # convert object to Lolliplot format
    lolliplotData <- toLolliplot(object, BSgenome=BSgenome, verbose=verbose)
    
    # retrieve the datasets from biomaRt
    mart <- retrieveMart(object=species, host=host, verbose=verbose)
    
    # annotate data with gene based on the ensembl transcript ID given, returns a list containing the mart
    lolliplotData <- annotateGene(lolliplotData, transcript=transcript, ensembl_mart=mart, verbose=verbose)
    
    # filter the data by gene to reduce the time for biomaRt queries
    lolliplotData <- filterByGene(lolliplotData, verbose=verbose)
    
    # annotate data with ensembl transcripts if necessary
    lolliplotData <- annotateTranscript(lolliplotData, ensembl_mart=mart, verbose=verbose)
    
    # grab the protien coordinates
    lolliplotData <- annotateProteinCoord(lolliplotData, ensembl_mart=mart, txdb=txdb, BSgenome=BSgenome, verbose=verbose)
    
    # filter data to only one transcript
    lolliplotData <- filterByTranscript(lolliplotData, transcript=transcript, verbose=verbose)
    
    # tabulate the frequency of observed mutations
    lolliplotData <- calcMutFreq(lolliplotData, verbose=verbose)
    
    # get the transcript size and domains
    proteinData <- constructTranscriptData(lolliplotData, ensembl_mart=mart, verbose=verbose)
      
    # set the mutations for the second tier
    lolliplotData <- setTierTwo(lolliplotData, proteinData=proteinData, emphasize=emphasize, verbose=verbose)
    
    # set the initial mutation heights (tier 1)
    lolliplotData <- setTierOne(lolliplotData, verbose=verbose)
    
    # set a label column
    lolliplotData <- addLabel(lolliplotData, verbose=verbose)
    
    # nest the protein domains
    proteinData <- setDomainHeights(proteinData, verbose=verbose)
    
    new("LolliplotData", primaryData=lolliplotData, geneData=proteinData)
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
LolliplotPlots <- function(object, labelAA, DomainPalette, MutationPalette, plotALayers, plotBLayers, verbose){
    
    # set up a density plot for the mutations
    PlotA <- buildDensityPlot(object, plotALayers=plotALayers, verbose=verbose)
    
    # set up the lolliplot
    PlotB <- buildLolliplot(object, labelAA=labelAA, DomainPalette=DomainPalette, MutationPalette=MutationPalette, plotBLayers=plotBLayers, verbose=verbose)
    
    new("LolliplotPlots", PlotA=PlotA, PlotB=PlotB)
}

################################################################################
########################### Accessor function definitions ######################

#' Helper function to get data from classes
#' 
#' @rdname getData-methods
#' @aliases getData
.getData_Lolliplot <- function(object, name=NULL, index=NULL, ...){
    
    if(is.null(name) & is.null(index)){
        memo <- paste("Both name and index are NULL, one must be specified!")
        stop(memo)
    }
    
    if(is.null(index)){
        index <- 0
    } else {
        if(index > 2){
            memo <- paste("index out of bounds")
            stop(memo)
        }
    }
    
    if(is.null(name)){
        name <- "noMatch"
    } else {
        slotAvailableName <- c("primaryData", "geneData")
        if(!(name %in% slotAvailableName)){
            memo <- paste("slot name not found, specify one of:", toString(slotAvailableName))
            stop(memo)
        }
    }
    
    if(name == "primaryData" | index == 1){
        data <- object@primaryData
    } else if(name == "geneData" | index == 2){
        data <- object@geneData
    }
    
    return(data)
}

#' @rdname getData-methods
#' @aliases getData
setMethod(f="getData",
          signature="LolliplotData",
          definition=.getData_Lolliplot)

#' @rdname getData-methods
#' @aliases getData
setMethod(f="getData",
          signature="Lolliplot",
          definition=.getData_Lolliplot)

#' Helper function to extract grobs from objects
#'
#' @rdname getGrob-methods
#' @aliases getGrob
#' @noRd
.getGrob_Lolliplot <- function(object, index=1, ...){
    if(index == 1){
        grob <- object@PlotA
    } else if(index == 2) {
        grob <- object@PlotB
    } else if(index == 3) {
        grob <- object@Grob
    } else {
        stop("Subscript out of bounds") 
    }
    return(grob)
}

#' @rdname getGrob-methods
#' @aliases getGrob
setMethod(f="getGrob",
          signature="LolliplotPlots",
          definition=.getGrob_Lolliplot)

#' @rdname getGrob-methods
#' @aliases getGrob
setMethod(f="getGrob",
          signature="Lolliplot",
          definition=.getGrob_Lolliplot)

#' @rdname drawPlot-methods
#' @aliases drawPlot
#' @importFrom grid grid.draw
#' @importFrom grid grid.newpage
#' @exportMethod drawPlot
setMethod(
    f="drawPlot",
    signature="Lolliplot",
    definition=function(object, ...){
        mainPlot <- getGrob(object, index=3)
        grid::grid.newpage()
        grid::grid.draw(mainPlot)
    }
)

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
          definition=function(object, verbose, ...){
              
              # get gene to filter by
              entrez <- object$biomaRt_entrez_gene_id[1]
              hgnc <- object$biomaRt_hgnc_symbol[1]
              ensembl <- object$biomaRt_ensembl_gene_id[1]
              
              # print status message
              if(verbose){
                  memo <- paste("Attempting to filter data by the following gene information: entrez =",
                                toString(entrez), ", hgnc =", toString(hgnc), ", ensembl =", toString(ensembl))
                  message(memo)
              }
              
              # check that the gene to filter is in the data
              if(!any(object$gene %in% c(entrez, hgnc, ensembl))){
                  memo <- paste("could not find the identifiers:",
                                toString(entrez), ", hgnc =", toString(hgnc), ", ensembl =", toString(ensembl),
                                "in the gene column of the data!")
                  stop(memo)
              }
              
              # find which identifier to use to filter
              if(any(object$gene %in% ensembl)){
                  filterGene <- ensembl
              } else if(any(object$gene %in% entrez)){
                  filterGene <- entrez
              } else {
                  filterGene <- hgnc
              }
              
              # another status message
              if(verbose){
                  memo <- paste("Choosing", toString(filterGene), "to filter data")
                  message(memo)
              }
              
              # perform the filtering and refactor the column
              index <- object$gene %in% filterGene
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
          definition=function(object, ensembl_mart, verbose, ...){
              
              # first check if both a transcript and proteinCoord column exist
              # we need both at this stage to ensure that the coordinates are accurate
              # and that we pull the correct domain data later on
              if(all(c("transcript", "proteinCoord") %in% colnames(object))){
                  
                  # next check to see if the ensembl transcript appears to be valid if so change the column name for downstream functions
                  if(any(grepl("^ENS\\D*T", object$transcript))){
                      
                      if(verbose){
                          memo <- paste("Found a transcript column with ensembl transcript annotations",
                                        "and amino acid coordinates... skipping transcript annotation")
                          message(memo)
                      }
                      object$ensembl_transcript_id <- object$transcript
                      object$transcript <- NULL
                      return(object)
                  }
              }
              
              # annotate the transcripts using the gene
              if(verbose){
                  memo <- paste("Annotating ensembl transcripts via biomaRt")
                  message(memo)
              }
             
              # Apply various filters using vector of values
              # biomaRt::listFilters(ensembl_mart)
              filters <- c("ensembl_gene_id")
              values <- c(as.character(object$biomaRt_ensembl_gene_id))
              
              # Select attributes to retrieve (protein domain, start, stop)
              # biomaRt::listAttributes(ensembl_mart)
              attributes <- c("ensembl_gene_id", "ensembl_transcript_id")
              
              # Retrieve data
              result <- biomaRt::getBM(attributes=attributes, filters=filters,
                                       values=values, mart=ensembl_mart)
              colnames(result) <- c("biomaRt_ensembl_gene_id", "biomaRt_ensembl_transcript_id")
              
              # add in the data and print a status message for what happened
              
              nrowOrig <- nrow(object)
              object <- merge(object, result, by=c("biomaRt_ensembl_gene_id"), allow.cartesian=TRUE)
              
              if(verbose){
                  memo <- paste("biomaRt found", toString(length(result$biomaRt_ensembl_transcript_id)),
                                "ensembl transcripts for gene id", toString(unique(result$biomaRt_ensembl_gene_id)),
                                ", expanding", nrowOrig, "entries to", nrow(object))
                  message(memo)
              }
              
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
          definition=function(object, transcript, ensembl_mart, verbose, ...){
             
              # perform quality checks on transcript
              if(is.null(transcript)){
                  memo <- paste("Parameter transcript is required, found NULL!")
                  stop(memo)
              }
              if(!is.character(transcript)){
                  memo <- paste("Parameter transcript is not of class character",
                                "found class:", toString(class(transcript)),
                                "attempting to coerce!")
                  warning(memo)
                  transcript <- as.character(transcript)
              }
              if(length(transcript) != 1){
                  transcript <- transcript[1]
                  memo <- paste("Parameter transcript has a length of", toString(length(transcript)),
                                "should be of length 1, using only the first element:", toString(transcript))
                  warning(memo)
              }

              # check that the value supplied to transcript appears to be an ensembl transcript
              if(grepl("^ENS\\D*G", transcript)){
                  memo <- paste("Parameter transcript appears to be an ensembl gene id,",
                                "should be an ensembl transcript id!")
                  stop(memo)
              }
              if(!grepl("^ENS\\D*T", transcript)){
                  memo <- paste("Parameter transcript does not appear to be an ensembl",
                                "transcipt ID, should follow the following form examples:",
                                "\"ENST0000\", \"ENSCCAT0000\", \"ENSMUST0000\", etc.")
                  stop(memo)
              }

              # Find the gene identifiers based on the transcript
              if(verbose){
                  memo <- paste("Retrieving gene identifiers based on the ensembl transcript id:", toString(transcript))
                  message(memo)
              }
              
              # Apply various filters using vector of values
              # biomaRt::listFilters(ensembl_mart)
              filters <- c("ensembl_transcript_id")
              values <- c(transcript)
              
              # Select attributes to retrieve (protein domain, start, stop)
              # biomaRt::listAttributes(ensembl_mart)
              attributes <- c("ensembl_gene_id", "hgnc_symbol", "entrezgene")
              
              # Retrieve data
              result <- biomaRt::getBM(attributes=attributes, filters=filters,
                                       values=values, mart=ensembl_mart)
              
              # add in the data and return the structure
              object$biomaRt_ensembl_gene_id <- result$ensembl_gene_id
              object$biomaRt_hgnc_symbol <- result$hgnc_symbol
              object$biomaRt_entrez_gene_id <- result$entrezgene
              
              return(object)
          })

#' @rdname annotateProteinCoord-methods
#' @aliases annotateProteinCoord
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @importFrom VariantAnnotation predictCoding
#' @importFrom Biostrings DNAStringSet
#' @importFrom biomaRt useMart
#' @importFrom biomaRt listDatasets
#' @importFrom biomaRt useDataset
#' @importFrom biomaRt getBM
#' @noRd
setMethod(f="annotateProteinCoord",
          signature="data.table",
          definition=function(object, ensembl_mart, txdb, BSgenome, verbose, ...){
              
              ##################################################################
              ############ if AA_change column is present ######################
              
              # first look for an amino acid coord column, if it's present get the
              # protein coordinates from there and do quality checks, we've already
              # checked for a valid ensembl transcript so we can be confident as to
              # the accuracy when pulling domains
              if(any(colnames(object) %in% "proteinCoord")){
                  if(all(is.numeric(object$proteinCoord))){
                   
                      return(object)
                  } else if(any(grepl("^p\\.", object$proteinCoord))){
                      
                      if(verbose){
                          countPNotation <- length(object$proteinCoord[grepl("^p\\.", object$proteinCoord)])
                          memo <- paste("Found", countPNotation,"entries with p. notation")
                          message(memo)
                      }
                      
                      if(!all(grepl("^p\\.", object$proteinCoord))){
                          
                          indexToRemove <- which(!grepl("^p\\.", object$proteinCoord))
                          
                          if(verbose){
                              memo <- paste("Removing", length(indexToRemove), "entries without p. notation, could these be intronic variants!")
                              message(memo)
                          }
                          object <- object[-indexToRemove,]
                  }
                  
                      object$proteinCoord <- as.numeric(gsub("p\\.[*a-zA-z]*(\\d+).*?$",
                                                             "\\1", object$proteinCoord,
                                                             perl=TRUE))
                      if(any(is.na(object$proteinCoord))){
                          invalidPos <- object[is.na(object$proteinCoord),c("chromosome", "start", "stop", "reference", "variant")]
                          invalidPos <- unique(do.call("paste", c(invalidPos, sep=":")))
                          memo <- paste("Failed to extract protein locations for the following coordinates:", toString(invalidPos),
                                        "please check p. notation for these positions!... removing these positions")
                          warning(memo)
                          object <- object[!is.na(object$proteinCoord),]
                      }
                      return(object)
                  } else {
                      if(verbose){
                          memo <- paste("Failed to extract protein coordinates",
                                        "attempting to annotate these using txdb object!")
                          message(memo)
                      }
                  }

              }
              
              ##################################################################
              ############ if AA_change column is not present ##################
              
              # check that required objects are present
              if(is.null(txdb) || is.null(BSgenome)){
                  memo <- paste("A valid txdb and BSgenome object are required when GenVisR must annotate",
                                "protein coordinates, but found null for these parameters.")
                  stop(memo)
              }
              
              # annotate the mutation coordinates given on the transcript level
              if(verbose){
                  memo <- paste("Attempting to annotate protein coordinates for transcripts")
              }

              # make sure BSgenome object matches txdb object
              seqBSgenome <- VariantAnnotation::seqlevels(BSgenome)
              seqtxdb <- VariantAnnotation::seqlevels(txdb)
              
              if(!all(seqBSgenome %in% seqtxdb) || !all(seqtxdb %in% seqBSgenome)){
                  mismatchSeqBSgenome <- seqBSgenome[!seqBSgenome %in% seqtxdb]
                  mismatchSeqtxdb <- seqtxdb[!seqtxdb %in% seqBSgenome]
                  mismatchSeq <- unique(c(mismatchSeqBSgenome, mismatchSeqtxdb))
                  
                  memo <- paste("The following entires are not in both the txdb and BSgenome objects supplied:",
                                toString(mismatchSeq), "This may affect protien annotations if mutations are in these regions!")
                  warning(memo)
              }
              
              # make sure lolliplotData matches the BSgenome and txdb object
              chrMismatch <- as.character(unique(object[!object$chromosome %in% seqtxdb,]$chromosome))
              
              if(length(chrMismatch) >= 1){
                  
                  memo <- paste("The following chromosomes do not match the supplied BSgenome object",
                                toString(chrMismatch))
                  warning(memo)
                  
                  # test if the chr mismatch is fixed by appending chr to chromosomes
                  chrMismatch_appendChr <- sum(as.character(unique(paste0("chr", object$chromosome))) %in% seqBSgenome)
                  if(chrMismatch_appendChr <= length(chrMismatch)){
                      memo <- paste("appending \"chr\" to chromosomes in attempt to fix mismatch with the BSgenome")
                      warning(memo)
                      object$chromosome <- paste0("chr", object$chromosome)
                  }
              }
              
              # convert data to a GRanges object
              gr1 <- GRanges(seqnames=object$chromosome, ranges=IRanges(start=object$start, end=object$stop), strand="*",
                             mcols=object[,c("biomaRt_ensembl_gene_id", "sample", "reference", "variant", "consequence", "biomaRt_ensembl_gene_id", "gene", "biomaRt_ensembl_transcript_id")])
              
              # get the protein coordinates
              gr1 <- VariantAnnotation::predictCoding(gr1, txdb, BSgenome, Biostrings::DNAStringSet(object$variant))
              gr1 <- data.table::as.data.table(gr1)
              
              # add in the TXNAME for the TXID from gr1
              txnameDT <- unique(data.table::as.data.table(select(txdb, gr1$TXID, "TXNAME", keytype="TXID")))
              txnameDT$TXID <- as.character(txnameDT$TXID)
              gr1$TXID <- as.character(gr1$TXID)
              object <- merge(gr1, txnameDT, by=c("TXID"))
              
              # for indels multiple amino acid positions can be reported, grab just the first one
              object$PROTEINLOC <- sapply(object$PROTEINLOC, function(x) x[1])
              
              # choose columns to keep
              keep <- c("mcols.sample", "seqnames", "start", "end", "mcols.reference", "mcols.variant",
                        "mcols.gene", "mcols.consequence", "mcols.biomaRt_ensembl_gene_id", "mcols.biomaRt_ensembl_transcript_id",
                        "PROTEINLOC", "TXNAME", "REFCODON", "VARCODON", "REFAA", "VARAA")
              object <- object[,keep,with=FALSE]
              colnames(object) <- c("sample", "chromosome", "start", "stop", "reference", "variant", "gene",
                                    "consequence", "ensembl_gene_id", "ensembl_transcript_id", "txdb.proteinCoord",
                                    "txdb.transcript", "txdb.refCodon", "txdb.varCodon", "txdb.refAA", "txdb.varAA")
              object$proteinCoord <- object$txdb.proteinCoord
              
              # Apply various filters using vector of values
              # biomaRt::listFilters(ensembl_mart)
              filters <- c("ucsc")
              values <- unique(c(object$txdb.transcript))
              
              # Select attributes to retrieve (protein domain, start, stop)
              # biomaRt::listAttributes(ensembl_mart)
              attributes <- c("ensembl_transcript_id", "ucsc")
              
              # Retrieve data
              result <- biomaRt::getBM(attributes=attributes, filters=filters,
                                       values=values, mart=ensembl_mart)
             
              if(nrow(result) == 0){
                  memo <- paste("An ensembl transcript id could not be found for the ucsc id:", toString(values),
                                "unable to continue... This has probably occurred because of a discrepancy between the txdb object",
                                "and the ensembl version queried.", "Please choose another transcript or alter the ensembl version queried (see the host parameter)!")
                  stop(memo)
              }
              
              if(length(unique(result$ensembl_transcript_id)) != length(unique(result$ucsc))){
                  memo <- paste("Retrieved protein coordinates are for the transcript", toString(unique(result$ucsc)),
                  "However GenVisR uses ensemble transcript id's", "there is not a 1 to 1 mapping between the two,",
                  "annotated coordinates may be inaccurate!")
                  warning(memo)
              }
              
              # use the biomaRt results to subset to only valid ensembl transcripts
              object$transcriptKey <- paste(object$ensembl_transcript_id, object$txdb.transcript, sep=":")
              result$transcriptKey <- paste(result$ensembl_transcript_id, result$ucsc, sep=":")
              object <- object[object$transcriptKey %in% result$transcriptKey,]
              object$transcriptKey <- NULL
              
              return(object)
          })

#' @rdname filterByTranscript-methods
#' @aliases filterByTranscript
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @noRd
setMethod(f="filterByTranscript",
          signature="data.table",
          definition=function(object, transcript, verbose, ...){
              
              # pick a transcript to filter on if none is given
              if(is.null(transcript)){
                  transcript <- as.character(unique(object$ensembl_transcript_id))[1]
                  if(verbose){
                      memo <- paste("parameter transcript is null, using transcript",
                                    toString(transcript), "available transcripts were:",
                                    toString(unique(object$ensembl_transcript_id)))
                      message(memo)
                  }
              } else {
                  # status message
                  if(verbose){
                      memo <- paste("Filtering data by transcript:", toString(transcript))
                      message(memo)
                  }
              }
              
              # perform quality checks
              if(!is.character(transcript)){
                  memo <- paste("input to transcript is not a character vector, attempting to coerce")
                  warning(memo)
              }
              if(length(transcript) > 1){
                  memo <- paste("input to transcript is not of length 1, using only the first element:", toString(transcript[1]))
                  warning(memo)
              }
              
              # check that the transcript is in the data
              if(!transcript %in% object$ensembl_transcript_id){
                  transcriptTmp <- as.character(unique(object$ensembl_transcript_id))[1]
                  memo <- paste("Did not find the specified transcript:", toString(transcript),
                                "in the data. Using:", toString(transcriptTmp), "instead!")
                  warning(memo)
                  transcript <- transcriptTmp
              }
              
              # subset the data.table
              object <- object[object$ensembl_transcript_id == transcript,]
              object$ensembl_transcript_id <- as.character(object$ensembl_transcript_id)
              
              return(object)
          })

#' @rdname calcMutFreq-methods
#' @aliases calcMutFreq
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @importFrom plyr count
#' @noRd
setMethod(f="calcMutFreq",
          signature="data.table",
          definition=function(object, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Calculating the frequency of mutations")
                  message(memo)
              }
              
              # make a temporary key making each mutation type unique, at this 
              # stage there should only be one transcript making genomic locations unique
              object$key <- paste(object$chromosome, object$start, object$stop, object$reference, object$variant, sep=":")
              
              # tabulate the frequencies
              mutationFreq <- plyr::count(object$key)
              colnames(mutationFreq) <- c("key", "mutationFreq")
              object <- merge(object, mutationFreq, by="key", all.x=TRUE)
              
              # remove the key column and return the object
              object$key <- NULL
              
              return(object)
          })

#' @rdname constructTranscriptData-methods
#' @aliases constructTranscriptData
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @importFrom plyr count
#' @noRd
setMethod(f="constructTranscriptData",
          signature="data.table",
          definition=function(object, ensembl_mart, verbose, ...){
              
              # make sure there's only one transcript
              if(length(unique(object$ensembl_transcript_id)) != 1){
                  memo <- paste("Number of transcripts does not equal one",
                                "unable to retrieve domains")
                  stop(memo)
              }
              
              # set the transcript
              transcript <- unique(object$ensembl_transcript_id)
              
              # print status message
              if(verbose){
                  memo <- paste("Retrieving protein domains for:", toString(transcript))
                  message(memo)
              }
              
              ############ Retrieve Domains ####################################
              
              # Apply various filters using vector of values
              #biomaRt::listFilters()
              filters <- c("ensembl_transcript_id")
              values <- as.list(c(transcript))
              
              # Select attributes to retrieve (protein domain, start, stop)
              #biomaRt::listAttributes()
              attributes <- c("interpro_description",
                              "interpro_short_description",
                              "interpro_start",
                              "interpro_end")
              
              # Retrieve data
              domain <- biomaRt::getBM(attributes=attributes, filters=filters,
                                       values=values, mart=ensembl_mart)
              domain <- as.data.table(domain)
              colnames(domain) <- c("interproDesc", "interproShortDesc", "proteinStart", "proteinStop")
              
              # check that domains were actually retrieved
              if(nrow(domain) < 1){
                  
                  if(verbose){
                      memo <- paste("biomaRt did not return protein domains")
                      message(memo)
                  }
                  
              } else {
                  domain$source <- "domain" 
              }
              
              
              ################ Retrieve transcript size ########################
              
              
              # Apply various filters using vector of values
              # biomaRt::listFilters()
              filters <- c("ensembl_transcript_id")
              values <- as.list(c(transcript))
              
              # Select attributes to retrieve (protein domain, start, stop)
              # biomaRt::listAttributes()
              attributes <- c("cds_length")
              
              # Retrieve data
              protein <- biomaRt::getBM(attributes=attributes, filters=filters,
                                        values=values, mart=ensembl_mart)
              protein <- as.data.table(protein)
              colnames(protein) <- "proteinStop"
              
              # make sure we have a protein length otherwise theres no gene to plot
              if(nrow(protein) < 1){
                  memo <- paste("biomaRt failed to return a CDS length for transcript:", toString(transcript))
                  warning(memo)
              } else {
                  protein$proteinStart <- 1
                  protein$source <- "protein"
                  protein$proteinStop <- protein$proteinStop/3
              }

              proteinData <- data.table::rbindlist(list(protein, domain), use.names=TRUE, fill=TRUE)
              
              return(proteinData)

          })

#' @rdname setTierOne-methods
#' @aliases setTierOne
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @noRd
setMethod(f="setTierOne",
          signature="data.table",
          definition=function(object, verbose, ...){
              
              # status message
              if(verbose){
                  memo <- paste("setting coordinates for first tier of mutations")
                  message(memo)
              }
              
              # split up by tier
              object <- split(object, by="tier")
              
              # set a base height for first tier
              object[["first"]]$height <- 1.2
              
              # make the height adjustment for first tier
              if(min(object[["first"]]$mutationFreq) != max(object[["first"]]$mutationFreq)){
                  # adjust the height of points based on their frequency as a scaling factor
                  a2 <- function(x, y, z){
                      ((1.4-1.2)*(x-y))/(z-y) + 1.2
                  }
                  object[["first"]]$height <- sapply(object[["first"]]$mutationFreq, a2, min(object[["first"]]$mutationFreq), max(object[["first"]]$mutationFreq))
              }
              
              # bind the object back together
              object <- rbindlist(object)
              
              # return the object with adjusted heights
              return(object)
          })

#' @rdname setTierTwo-methods
#' @aliases setTierTwo
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @importFrom plyr count
#' @importFrom plyr mapvalues
#' @noRd
setMethod(f="setTierTwo",
          signature="data.table",
          definition=function(object, proteinData, emphasize, verbose, ...){
              
              # status message
              if(verbose){
                  memo <- paste("setting coordinates for second tier of mutations")
                  message(memo)
              }
              
              #################### set level 2 mutation base height ############
              
              # by default grab the top 10% of most frequently mutated events for the second tier
              if(!is.null(emphasize)){
                  
                  # perform quality checks on emphasize
                  if(!is.numeric(emphasize)){
                      memo <- paste("values in emphasize are not integers, attempting to coerce.")
                      warning(memo)
                      emphasize <- as.integer(emphasize)
                  }
                  notFound <- emphasize[!emphasize %in% object$proteinCoord]
                  if(length(notFound) > 0){
                      memo <- paste("The following coordinates given in emphasize:", toString(notFound),
                                    "were not found!")
                      warning(memo)
                  }
                  
                  # set everything in emphasize to the second tier
                  object$tier <- ifelse(object$proteinCoord %in% emphasize, "second", "first")
                  
              } else {
                  
                  # get all the protein coordinates
                  proteinLoc <- plyr::count(object$proteinCoord)
                  colnames(proteinLoc) <- c("Coord", "Frequency")
                  proteinLoc <- proteinLoc[order(proteinLoc$Frequency, proteinLoc$Coord),]
                  
                  # find how many are in the to ten percent
                  top <- ceiling(.1 * nrow(proteinLoc))
                  topProteinLoc <- tail(proteinLoc$Coord, top)
                  
                  # set the tier
                  object$tier <- ifelse(object$proteinCoord %in% topProteinLoc, "second", "first")
              }
              
              # set base height
              object[object$tier == "second","height"] <- 2.6
              
              ############### set level 2 A mutation heights ###################
              
              # split up object by the protein coordinate
              object <- split(object, by="proteinCoord")
              
              # set up anonymous function to set level 2A heights
              a <- function(x){
                  
                  # if a coordinate is in the first tier of mutations do nothing
                  if(x$tier[1] == "first"){
                      return(x)
                  }
                  
                  # order by the frequency of mutation events for the coordinate
                  x$key <- paste(x$chromosome, x$start, x$stop, x$reference, x$variant, sep=":")
                  x <- x[order(-x$mutationFreq),]
                  
                  # consecutively number each unique key
                  x <- x[, tmp := .GRP, by=list(key)][]
                  
                  # convert each consequitive number to a scaling factor i.e 1 -> 0, 2 -> .1 etc.
                  # then use those scaled as a height adjustment for level 2A
                  x$tmp <- x$tmp/10 - .1
                  x$height <- x$tmp + x$height
                  x$tmp <- NULL
                  x$key <- NULL
                  
                  return(x)
              }
              
              object <- lapply(object, a)
              object <- rbindlist(object)
              
              ####### adjust level 2 mutations to span across the gene #########
              
              # find values for those points in second tier based on the protein
              protein <- proteinData[proteinData$source == "protein",]
              from <- protein$proteinStart + protein$proteinStop/top/2
              to <- protein$proteinStop
              by <- protein$proteinStop/top
              newValues <- seq(from=from, to=to, by=by)
              
              # map the adjusted coordinate values
              object$proteinCoordAdj <- object$proteinCoord
              object <- object[order(object$proteinCoordAdj),]
              object$proteinCoordAdj <- plyr::mapvalues(object$proteinCoordAdj, unique(object[object$tier == "second",]$proteinCoordAdj), newValues)
              
              # return the object with adjusted heights
              return(object)
          })

#' @rdname toLolliplot-methods
#' @aliases toLolliplot
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @noRd
setMethod(f="toLolliplot",
          signature="data.table",
          definition=function(object, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("converting", class(object)[1], "to expected Lolliplot format")
                  message(memo)
              }
              
              # check for the right columns
              expecCol <- c("chromosome", "start", "stop", "reference", "variant", "sample", "gene")
              if(!all(expecCol %in% colnames(object))){
                  missingCol <- expecCol[!expecCol %in% colnames(object)]
                  memo <- paste("The following required columns were missing from the input:",
                                toString(missingCol))
                  stop(memo)
              }
              
              # return the object
              return(object)
          })

#' @rdname toLolliplot-methods
#' @aliases toLolliplot
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @noRd
#' @importFrom data.table as.data.table
setMethod(f="toLolliplot",
          signature="data.frame",
          definition=function(object, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("converting", class(object), "to data.table")
                  message(memo)
              }
              
              # first convert to data.table
              object <- as.data.table(object)
              
              # then run toLolliplot on the new data.table
              object <- toLolliplot(object, verbose=verbose)
              
              return(object)
          })

#' @rdname retrieveMart-methods
#' @aliases retrieveMart
#' @param object Object of class character giving a species
#' @param verbose Boolean specifying if status messages should be reported
#' @noRd
setMethod(f="retrieveMart",
          signature="character",
          definition=function(object, host, verbose, ...){
              
              # status message
              if(verbose){
                  memo <- paste("attempting to retrieve biomaRt dataset for",
                                "species:", toString(object))
                  message(memo)
              }
              
              # set the species
              species <- object
              
              # load the appropriate mart
              ensembl_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                                               host=host)
              
              # select proper data set given regexp print warnings if unexpected out occur
              dataset <- biomaRt::listDatasets(ensembl_mart)$dataset
              index <- which(grepl(species, dataset))
              if(length(index)>1)
              {
                  memo <- paste(toString(species), " Matches more than one dataset for the",
                                " ensembl mart, please specify the full dataset from this list:",
                                toString(dataset))
                  stop(memo)
              } else if(length(index)==0) {
                  memo <- paste(toString(species), " does not appear to be supported by biomaRt",
                                "available datasets are:", toString(dataset))
                  stop(memo)
              }
              
              ensembl_mart <- biomaRt::useDataset(as.character(dataset[index]),
                                                  mart=ensembl_mart)
              
              return(ensembl_mart)
          })

#' @rdname buildDensityPlot-methods
#' @aliases Lolliplot
#' @param object Object of class LolliplotData
#' @return gtable object containg a density plot for the lolliplot
#' @noRd
#' @import ggplot2
setMethod(f="buildDensityPlot",
          signature="LolliplotData",
          definition=function(object, plotALayers, verbose=verbose){
              
              # extract the data we need
              primaryData <- getData(object, name="primaryData")
              proteinData <- getData(object, name="geneData")
              
              # make a pseudo data set containing the min/max coordinates possible
              min1 <- min(c(proteinData$proteinStart, primaryData$proteinCoordAdj, primaryData$proteinCoord))
              max1 <- max(c(proteinData$proteinStop, primaryData$proteinCoordAdj, primaryData$proteinCoord))
              coordData <- as.data.frame(c(min1, max1)) 
              colnames(coordData) <- c("start")
              
              # print status message
              if(verbose){
                  memo <- paste("Building the Density plot")
                  message(memo)
              }
              
              # perform quality checks
              if(!is.null(plotALayers)){
                  if(!is.list(plotALayers)){
                      memo <- paste("plotALayers is not a list")
                      stop(memo)
                  }
                  
                  if(any(!unlist(lapply(plotALayers, function(x) ggplot2::is.ggproto(x) | ggplot2::is.theme(x) | is(x, "labels"))))){
                      memo <- paste("plotALayers is not a list of ggproto or ",
                                    "theme objects... setting plotALayers to NULL")
                      warning(memo)
                      plotALayers <- NULL
                  }
              }
              
              ############# start building the plot ############################
              
              # base theme
              baseTheme <- theme_bw()
              
              # plot theme
              plotTheme <- theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                                 panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
                                 axis.title.x=element_blank())
              
              # define coord
              boundaries <- geom_blank(data=coordData, aes_string(x='start'))
              
              # plot labels
              
              # define the geom
              plotGeom <- geom_density(mapping=aes_string(x="proteinCoord"), fill="lightcyan4", color="lightcyan4", adjust=.75)
              
              # define the plot
              densityPlot <- ggplot(data=primaryData)
              
              # combine everything
              densityPlot <- densityPlot + plotGeom + baseTheme + plotTheme + boundaries + plotALayers 
              
              # convert to grob
              densityGrob <- ggplotGrob(densityPlot)
              return(densityGrob)
          })

#' @rdname buildLolliplot-methods
#' @aliases Lolliplot
#' @param object Object of class LolliplotData
#' @return gtable object containg a lolliplot
#' @noRd
#' @import ggplot2
setMethod(f="buildLolliplot",
          signature="LolliplotData",
          definition=function(object, labelAA, plotBLayers, DomainPalette, MutationPalette, verbose=verbose){
              
              # extract the data we need
              primaryData <- getData(object, name="primaryData")
              proteinData <- getData(object, name="geneData")
              
              # split the data up for convenience
              domainData <- proteinData[proteinData$source == "domain",]
              proteinData <- proteinData[proteinData$source == "protein",]
              firstTier <- primaryData[primaryData$tier == "first",]
              secondTier <- primaryData[primaryData$tier == "second",]
              
              # print status message
              if(verbose){
                  memo <- paste("Building the lolliplot")
                  message(memo)
              }
              
              # perform quality checks
              if(!is.null(plotBLayers)){
                  if(!is.list(plotBLayers)){
                      memo <- paste("plotBLayers is not a list")
                      stop(memo)
                  }
                  
                  if(any(!unlist(lapply(plotBLayers, function(x) ggplot2::is.ggproto(x) | ggplot2::is.theme(x) | is(x, "labels"))))){
                      memo <- paste("plotBLayers is not a list of ggproto or ",
                                    "theme objects... setting plotBLayers to NULL")
                      warning(memo)
                      plotBLayers <- NULL
                  }
              }
              
              # set plot palette for domains and assign random colors if there are not enough for each domain
              if(is.null(DomainPalette)){
                  plotPaletteDomain <- c("#94C873", "#73C8AC", "#73C3C8", "#73A0C8", "#7379C8",
                                         "#A573C8", "#C8737F", "#D8D444", "#94AEA1", "#547C81",
                                         "#515587", "#875151", "#C3923D")
              } else {
                  plotPaletteDomain <- DomainPalette
              }
              
              if(length(plotPaletteDomain) < length(unique(domainData$interproDesc))){
                  
                  # if user supplied palette with not enough colors give warning
                  if(!is.null(DomainPalette)){
                      memo <- paste("Found", length(plotPaletteDomain), "values in palette but there are", length(unique(domainData$interproDesc)),
                                    "protein domains... adding random colors to complete Domainpalette!")
                      warning(memo)
                  }

                  newCol <- grDevices::colors(distinct=TRUE)[!grepl("^gray", grDevices::colors(distinct=TRUE))]
                  num2Choose <- length(unique(domainData$interproDesc) )- length(plotPaletteDomain)
                  newCol <- sample(newCol, num2Choose)
                  plotPaletteDomain <- c(plotPaletteDomain, newCol)
              }
              
              # set plot palette for mutations and assign random colors if there are not enough for each mutationt type
              MutationPalette
              if(is.null(MutationPalette)){
                  plotPaletteMutation <- c("#0c457d", "#87ceff", "#ffcc5c", "#e8702a", "#9a8262", "#725c91", "#7e9d97", "#9f4f4f", "#4f6333")
              } else {
                  plotPaletteMutation <- MutationPalette
              }
              
              if(length(plotPaletteMutation) < length(unique(primaryData$consequence))){
                  
                  # if user supplied palette with not enough colors give warning
                  if(!is.null(MutationPalette)){
                      memo <- paste("Found", length(plotPaletteMutation), "values in palette but there are", length(unique(primaryData$consequence)),
                                    "mutation types... adding random colors to complete Mutationpalette!")
                      warning(memo)
                  }
                  
                  newCol <- grDevices::colors(distinct=TRUE)[!grepl("^gray", grDevices::colors(distinct=TRUE))]
                  num2Choose <- length(unique(primaryData$consequence) )- length(plotPaletteMutation)
                  newCol <- sample(newCol, num2Choose)
                  plotPaletteMutation <- c(plotPaletteMutation, newCol)
              }
              
              ############# start building the plot ############################
              
              # set up the protein plot
              protein <- geom_rect(data=proteinData, mapping=aes_string(xmin="proteinStart", xmax="proteinStop", ymin="minHeight", ymax="maxHeight"))
              
              # set up the domains
              domain <- geom_rect(data=domainData, mapping=aes_string(xmin="proteinStart", xmax="proteinStop", ymin="minHeight", ymax="maxHeight", fill="interproDesc"))
              
              # set up the lolli segments
              segment1 <- geom_segment(data=secondTier, mapping=aes_string(x="proteinCoordAdj", y="height", xend="proteinCoordAdj", yend=2.2), color="red")
              segment2 <- geom_segment(data=secondTier, mapping=aes_string(x="proteinCoordAdj", y=2.2, xend="proteinCoord", yend=1.5), color="red")
              segment3 <- geom_segment(data=secondTier, mapping=aes_string(x="proteinCoord", y=1.5, xend="proteinCoord", yend=1), color="red")
              segment4 <- geom_segment(data=firstTier, mapping=aes_string(x="proteinCoord", y="height", xend="proteinCoord", yend=1))
              
              # set up the text
              if(labelAA){
                  nudgeFactor <- .01 * proteinData$proteinStop[1]
                  plotText <- geom_text(data=secondTier, mapping=aes_string(x="proteinCoordAdj", y="height", label="label"), position=position_nudge(x=nudgeFactor), hjust=0)
              } else {
                  plotText <- geom_blank()
              }
              
              # set up the points
              points <- geom_point(data=primaryData, aes_string(x="proteinCoordAdj", y="height", size="mutationFreq", color="consequence"), shape=19)
              
              # set base theme
              baseTheme <- theme_bw()
              
              # set plot theme
              plotTheme <- theme(legend.position="bottom", axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                                 axis.title.y=element_blank(), panel.grid.major.y=element_blank(),
                                 panel.grid.minor.y=element_blank())
              
              # set colors for plot
              plotDomainPalette <- scale_fill_manual(values=plotPaletteDomain)
              plotMutationPalette <- scale_color_manual(values=plotPaletteMutation)
              
              # set plot labels
              plotLabel <- labs(x="Amino Acid Position", fill="Interpro\nProtein Domain", size="Mutation\nFrequency", color="Mutation\nConsequence")
              
              # define legend guide
              legGuide <- guides(fill=guide_legend(label.position="right", title.position="top"),
                                 color=guide_legend(label.position="right", title.position="top", ncol=3),
                                 size=guide_legend(label.position="right", title.position="top", ncol=2))
              
              # define the plot
              lolliplotPlot <- ggplot()
              
              # combine everything
              lolliplotPlot <- lolliplotPlot + protein + domain + segment1 + segment2 + segment3 + segment4 + points +
                               baseTheme + plotTheme + plotLabel + plotDomainPalette + plotMutationPalette + plotText +
                               legGuide + scale_size(range=c(5, 10))

              # convert to grob
              lolliplotGrob <- ggplotGrob(lolliplotPlot)
              
              return(lolliplotGrob)
          })

#' @rdname arrangeLolliplotPlot-methods
#' @aliases arrangeLolliplotPlot
#' @param object Object of class LolliplotData
#' @return gtable object containg a lolliplot
#' @noRd
#' @importFrom gridExtra arrangeGrob
setMethod(f="arrangeLolliplotPlot",
          signature="LolliplotPlots",
          definition=function(object, sectionHeights, verbose=verbose){
              
              # grab the data we need
              plotA <- getGrob(object, 1)
              plotB <- getGrob(object, 2)
              
              # Obtain the max width for relevant plots
              plotList <- list(plotA, plotB)
              plotList <- plotList[lapply(plotList, length) > 0]
              plotWidths <- lapply(plotList, function(x) x$widths)
              maxWidth <- do.call(grid::unit.pmax, plotWidths)
              
              # Set the widths for all plots
              for(i in 1:length(plotList)){
                  plotList[[i]]$widths <- maxWidth
              }
              
              # set section heights based upon the number of sections
              defaultPlotHeights <- c(.25, .75)
              
              if(is.null(sectionHeights) & length(plotList) == length(defaultPlotHeights)){
                  sectionHeights <- defaultPlotHeights
              } else if(length(sectionHeights) != length(plotList)){
                  memo <- paste("There are", length(sectionHeights), "section heights provided",
                                "but", length(plotList), "vertical sections...",
                                "using default values!")
                  warning(memo)
                  sectionHeights <- defaultPlotHeights
              } else if(!all(is.numeric(sectionHeights))) {
                  memo <- paste("sectionHeights must be numeric... Using",
                                "default values!")
                  warning(memo)
                  sectionHeights <- defaultPlotHeights
              }
              
              # arrange the final plot
              finalPlot <- do.call(gridExtra::arrangeGrob, c(plotList, list(ncol=1, heights=sectionHeights)))
              
              return(finalPlot)
          })

#' @rdname setDomainHeights-methods
#' @aliases setDomainHeights
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @noRd
setMethod(f="setDomainHeights",
          signature="data.table",
          definition=function(object, verbose, ...){
              
              # status message
              if(verbose){
                  memo <- paste("Nesting", toString(nrow(object) - 1), "protein domains")
                  message(memo)
              }
              
              # add in min and max height columns for the protein and it's domains
              object$minHeight <- -1
              object$maxHeight <- 1
              
              # separate out the protein and the domains
              protein <- object[object$source == "protein",]
              domain <- object[object$source == "domain",]
              
              # attempt to simplify domain coordinates if any are overlapping
              domain <- split(domain, by="interproShortDesc")
              a <- function(x){
                  
                  # we can convert to granges ad use reduce to do this for us
                  gr1 <- GRanges(seqnames=1, ranges=IRanges(start=x$proteinStart, end=x$proteinStop), strand="*")
                  gr1 <- reduce(gr1)
                  gr1 <- as.data.table(gr1)
                  keep <- c("start", "end")
                  gr1 <- gr1[, keep, with=FALSE]
                  colnames(gr1) <- c("proteinStart", "proteinStop")
                  
                  # reconstruct the data frame
                  x <- unique(x[,-c("proteinStop", "proteinStart")])
                  x <- cbind(gr1, x)
                  
                  return(x)
              }
              domain <- lapply(domain, a)
              domain <- rbindlist(domain)
              
              # determine which regions are overlapping and annotate which nest domain is
              # sort on start
              domain$proteinStart <- as.numeric(domain$proteinStart)
              domain$proteinStop <- as.numeric(domain$proteinStop)
              domain <- domain[order(domain$proteinStart),]
              
              # annotate nests
              nest <- rep(1, nrow(domain))
              for(i in 2:nrow(domain)) {
                  if(domain$proteinStart[i] < domain$proteinStop[i-1]) {
                      # if start occurs before last stop
                      nest[i] <- nest[i-1]+1
                  } else {
                      # if start occurs after last stop
                      # get number of domains underneath
                      layer <- domain$proteinStart[i] < domain$proteinStop[seq(1,i-1)]
                      
                      # take care of edge case where no layer is true
                      if(all(!layer)){
                          nest[i] <- 1
                      } else {
                          nest[i] <- max(nest[which(layer==TRUE)])+1
                      }
                      
                  }
              }
              
              # add this nest information to the data frame
              domain$nest <- nest
              
              # combine gene and domain information
              protein$nest <- 1
              protein <- rbind(protein, domain)
              
              # annotate display heights based on nesting and make sure coord are numeric
              protein$minHeight <- 1/(as.numeric(protein$nest))
              protein$maxHeight <- -1/(as.numeric(protein$nest))
              
              return(protein)
          })

#' @rdname addLabel-methods
#' @aliases addLabel
#' @param object Object of class data.table
#' @param verbose Boolean specifying if status messages should be reported
#' @noRd
setMethod(f="addLabel",
          signature="data.table",
          definition=function(object, verbose, ...){
              
              # status message
              if(verbose){
                  memo <- paste("Setting labels for mutations")
                  message(memo)
              }
              
              # make a label column if there isn't one already
              if(!"label" %in% colnames(object)){
                  
                  # split into insertions, deletion, snvs and add a label
                  insIndex <- which(grepl("[-0]", object$reference))
                  ins <- object[insIndex,]
                  ins$label <- paste0("p.", ins$proteinCoord, "Ins")
                  
                  delIndex <- which(grepl("[-0]", object$variant))
                  del <- object[delIndex,]
                  del$label <- paste0("p.", del$proteinCoord, "Del")
                  
                  snvsIndex <- which(grepl("[ACGTacgt]", object$reference) & grepl("[ACGTacgt]", object$variant))
                  snvs <- object[snvsIndex,]
                  snvs$label <- paste0("p.", snvs$txdb.refAA, snvs$proteinCoord, snvs$txdb.varAA)
                  
                  other <- object[!rownames(object) %in% c(snvsIndex, delIndex, insIndex),]
                  other$label <- NA
                  
                  object <- data.table::rbindlist(list(snvs, ins, del, other))
                  
              }
              
              # reformat the label column so there is only one label per x,y coordinate
              dupLabIndex <- which(duplicated(object[,c("height", "proteinCoordAdj", "label")]))
              object[dupLabIndex, c("label")] <- NA
              
              return(object)
              
          })

