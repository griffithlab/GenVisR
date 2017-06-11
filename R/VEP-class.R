#' Class VEP
#' 
#' An S4 class for Variant Effect Predictor input.
#' @name VEP-class
#' @rdname VEP-class
#' @slot path Character string specifying the paths of the VEP files read in.
#' @slot version Numeric value specifying the version of VEP used.
#' @slot vepObject vep object which inherits from VEP_Virtual class.
#' @exportClass VEP
#' @include VEP_Virtual-class.R
#' @import methods
setClass("VEP",
         representation=representation(path="character",
                                       version="numeric",
                                       vepObject="VEP_Virtual"))

#' Initalizer method for the VEP container class
#' 
#' @name VEP
#' @rdname VEP-class
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
setMethod(
    f="initialize",
    signature="VEP",
    definition=function(.Object, path, data, version, verbose){
        cat("!!!!! VEP~Initalizer !!!!!\n")
        
        # Grab all files and assign to slot
        vepFiles <- Sys.glob(path)
        .Object@path <- vepFiles
        
        # anonymous function to read in files
        a <- function(x, verbose){
            # detect OS and remove slashes and extension
            if(.Platform$OS.type == "windows"){
                sampleName <- gsub("(.*/)||(.*\\\\)", "", x)
                sampleName <- gsub("\\.[^.]+$", "", x)
            } else {
                sampleName <- gsub("(.*/)", "", x)
                sampleName <- gsub("\\.[^.]+$", "", sampleName)
            }
            # read the header
            header <- readLines(con=x, n=400)
            header <- header[grepl("^##", header)]
            # find where headers stop and read the data
            skip <- length(header)
            vepData <- suppressWarnings(data.table::fread(input=x,
                                                          stringsAsFactors=TRUE,
                                                          verbose=verbose,
                                                          skip=skip))
            # set sample if it's not already in the data table
            if(any(colnames(vepData) %in% "sample")){
                return(vepData)
            } else {
                vepData$sample <- sampleName
                return(list("data"=vepData, "header"=header))
            }
        }
       
        # aggregate data into a single data table if necessary
        if(length(vepFiles) == 0){
            memo <- paste("No files found using:", path)
            stop(memo)
        } else {
            # Read in the information
            vepInfo <- lapply(vepFiles, a, verbose)
            
            # extract header and data information
            vepHeader <- lapply(vepInfo, function(x) x[["header"]])
            vepData <- lapply(vepInfo, function(x) x[["data"]])
            
            # aggregate the data
            vepData <- data.table::rbindlist(vepData, fill=TRUE)
        } 
        
        # grab the version and assign it
        a <- function(x){
            # find the element which defines the VEP version
            x <- x[grepl("VARIANT EFFECT PREDICTOR", x)]
            
            # extract the version
            x <- regmatches(x,regexpr("[0-9]+\\.*[0-9]*",x))
            
            if(length(x) != 1){
                memo <- paste("Expected 1 entry for VEP version, found:",
                              length(x), "using", as.numeric(x[1]))
                warning(memo)
            }
            return(as.numeric(x[1]))
        }
        if(version == "auto"){
            version <- lapply(vepHeader, a)
            version <- unique(unlist(version))
            if(length(version) > 1){
                version <- version[1]
                memo <- paste("Expect 1 version, the following versions were",
                              "found:", toString(version), "Using version",
                              version, "for parsing!")
                warning(memo)
            } else if(length(version) == 0){
                memo <- paste("Cannot infer version from vep headers",
                              "no versions found!")
                stop(memo)
            }
            .Object@version <- version
        } else {
            .Object@version <- version
        }
        
        # assign the vepData to it's slot
        if(version >= 88 || version < 89){
            .Object@vepObject <- VEP_v88(vepData=vepData, vepHeader=vepHeader)
        } else {
            memo <- paste("Currently only VEP version 88 is supported, make a",
                          "feature request on",
                          "https://github.com/griffithlab/GenVisR!")
            stop(memo)
        }
        
        return(.Object)
    })

#' Constructor for the VEP container class.
#' 
#' @name VEP
#' @rdname VEP-class
#' @param path String specifying the path to a VEP annotation file. Can accept
#' wildcards if multiple VEP annotation files exist (see details).
#' @param version String specifying the version of the VEP files, Defaults to
#' auto which will look for the version in the header.
#' @param verbose Boolean specifying if progress should be reported while
#' reading in the VEP files.
#' @details When specifying a path to a VEP annotation file the option exist to
#' either specify the full path to an annotation file or to use wildcards to
#' specify multiple files. When specifying a full path the initalizer will check
#' if a column named "sample" containg the relevant sample for each row exists.
#' If such a column is not found the initalizer will assume this file
#' corresponds to only one sample and populate a sample column accordingly.
#' Alternatively if multiple files are specified at once using a wildcard, the
#' initalizer will aggregate all the files and use the file names minus any
#' extension to populate sample names.
#' @export
VEP <- function(path, version="auto", verbose=FALSE){
    cat("!!!!! VEP~Constructor !!!!!\n")
    new("VEP", path=path, version=version, verbose=verbose)}

#' @rdname writeData-methods
#' @aliases writeData,VEP
#' @param object Object of class VEP.
#' @param file Character string specifying a file to send output to.
#' @param sep Delimiter used when writing output, defaults to tab.
setMethod(f="writeData",
          signature="VEP",
          definition=function(object, file, ...){
              writeData(object@vepObject, file, sep="\t")
          })

#' @rdname setMutationHierarchy-methods
#' @aliases setMutationHierarchy,VEP
#' @noRd
#' @importFrom data.table data.table
#' @importFrom data.table setDT
#' @importFrom data.table rbindlist
setMethod(f="setMutationHierarchy",
          signature="VEP",
          definition=function(object, mutationHierarchy, verbose, ...){
              # set the mutation hierarchy if a custom hierarchy is unspecified
              if(is.null(mutationHierarchy)){
                  mutationHierarchy$mutation <- c("transcript_ablation", "splice_acceptor_variant",
                                                  "splice_donor_variant", "stop_gained",
                                                  "frameshift_variant", "stop_lost", "start_lost",
                                                  "transcript_amplification", "inframe_insertion",
                                                  "inframe_deletion", "missense_variant",
                                                  "protein_altering_variant", "splice_region_variant",
                                                  "incomplete_terminal_codon_variant", "stop_retained_variant",
                                                  "synonymous_variant", "coding_sequence_variant",
                                                  "mature_miRNA_variant", "5_prime_UTR_variant",
                                                  "3_prime_UTR_variant", "non_coding_transcript_exon_variant",
                                                  "intron_variant", "NMD_transcript_variant",
                                                  "non_coding_transcript_variant", "upstream_gene_variant",
                                                  "downstream_gene_variant", "TFBS_ablation",
                                                  "TFBS_amplification", "TF_binding_site_variant",
                                                  "regulatory_region_ablation", "regulatory_region_amplification",
                                                  "feature_elongation", "regulatory_region_variant",
                                                  "feature_truncation", "intergenic_variant")
                  
                  mutationHierarchy$color <- c("#bd3d3d", "#d45a73", "#cd3300",
                                               "#d24719", "#895878", "#e18466",
                                               "#bd92cc", "#9e55a7", "#7d4281",
                                               "#a6cea9", "#78b3d2", "#00a0b0",
                                               "#4f372d", "#c6c386", "#739475",
                                               "#ffdb86", "#ddb554", "#d69f35",
                                               "#f98107", "#383838", "#5acb7f",
                                               "#1e81b5", "#bcd9e8", "#ab763d",
                                               "#a72ca2", "#f4a460", "#0000ff",
                                               "#00510a", "#2a7694", "#00c9bb",
                                               "#8c860e", "#687233", "#6f3333",
                                               "#626161", "#db9294")
                  mutationHierarchy <- data.table::data.table("mutation"=mutationHierarchy$mutation,
                                                              "color"=mutationHierarchy$color)
              }
              
              # check that mutationHiearchy is a data table
              if(!any(class(mutationHierarchy) %in% "data.table")){
                  memo <- paste("mutationHiearchy is not an object of class",
                                "data.table, attempting to coerce.")
                  warning(memo)
                  mutationHierarchy <- data.table::setDT(mutationHierarchy)
              }
              
              # check for the correct columns
              correctCol <- c("mutation", "color")
              if(!all(correctCol %in% colnames(mutationHierarchy))){
                  missingCol <- correctCol[!correctCol %in% colnames(mutationHierarchy)]
                  memo <- paste("The correct columns were not found in",
                                "mutationHierarchy, missing", toString(missingCol))
                  stop(memo)
              }
              
              # check that all mutations are specified
              if(!all(object@vepObject@mutation$trv_type %in% mutationHierarchy$mutation)){
                  missingMutations <- unique(object@vepObject@mutation$trv_type[!object@vepObject@mutation$trv_type %in% mutationHierarchy$mutation])
                  memo <- paste("The following mutations were found in the",
                                "input however were not specified in the",
                                "mutationHierarchy!", toString(missingMutations),
                                "adding these in as least important and",
                                "assigning random colors!")
                  warning(memo)
                  newCol <- colors(distinct=TRUE)[!grepl("^gray", colors(distinct=TRUE))]
                  tmp <- data.table::data.table("mutation"=missingMutations,
                                                "color"=sample(newCol, length(missingMutations)))
                  mutationHierarchy <- data.table::rbindlist(list(mutationHierarchy, tmp), use.names=TRUE, fill=TRUE)
              }
              
              # add in a pretty print mutation labels
              mutationHierarchy$label <- gsub("_", " ", mutationHierarchy$mutation)
              mutationHierarchy$label <-  gsub("'", "' ", mutationHierarchy$mutation)
              
              # check for duplicate mutations
              if(any(duplicated(mutationHierarchy$mutation))){
                  duplicateMut <- mutationHierarchy[duplicated(mutationHierarchy$mutation),"mutation"]
                  memo <- paste("The mutation type",toString(duplicateMut),
                                "was duplicated in the supplied mutationHierarchy!")
                  mutationHierarchy <- mutationHierarchy[!duplicated(mutationHierarchy$mutation),]
              }
              
              # print status message
              if(verbose){
                  memo <- paste("Setting the hierarchy of mutations from most",
                                "to least deleterious and mapping to colors:",
                                toString(mutationHierarchy$mutation))
                  message(memo)
              }
              
              return(mutationHierarchy)
          })

#' @rdname toWaterfall-methods
#' @aliases toWaterfall,VEP
#' @noRd
setMethod(f="toWaterfall",
          signature="VEP",
          definition=function(object, hierarchy, labelColumn, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Converting", class(object),
                                "to expected waterfall format")
                  message(memo)
              }
              
              # grab the mutation hierarchy
              hierarchy <- hierarchy@MutationHierarchy
              
              # grab the sample, mutation, gene columns and set a label
              sample <- object@vepObject@sample
              mutation <- object@vepObject@mutation[,"Consequence"]
              gene <- object@vepObject@meta[,"SYMBOL"]
              label <- NA
              
              # if a label column exists and is proper overwrite the label variable
              if(!is.null(labelColumn)){
                  if(length(labelColumn) != 1) {
                      memo <- paste("Parameter \"labelColumn\" must be of length 1!",
                                    "Found length to be", length(labelColumn))
                      warning(memo)
                      next
                  } else if(labelColumn %in% colnames(object@gmsObject@meta)){
                      memo <- paste("Did not find column:", labelColumn,
                                    "in the meta slot of the vepObject! Valid",
                                    "names are:", colnames(getMeta(object)))
                      warning(memo)
                      next
                  } else {
                      label <- object@vepObject@meta[,labelColumn]
                  }
              }
              
              # combine all columns into a consistent format
              waterfallFormat <- cbind(sample, gene, mutation, label)
              colnames(waterfallFormat) <- c("sample", "gene", "mutation", "label")
              
              # make a temporary ID column for genomic features to collapse on
              # this will ensure the mutation burden/frequency plot will be accurate
              waterfallFormat$key <- paste0(object@vepObject@position$Location, ":",
                                            object@vepObject@mutation$Allele, ":",
                                            object@vepObject@sample$sample)
              rowCountOrig <- nrow(waterfallFormat)
              
              # cast data into form where mutation column is not comma delimited
              # the hierarchy will sort out any duplicates
              waterfallFormat <- waterfallFormat[, strsplit(as.character(mutation), ",", fixed=TRUE),
                                                 by = .(sample, gene, label, mutation, key)][,.(Consequence = V1, sample, gene, label, key)]
              
              # order the data based on the mutation hierarchy,
              # remove all duplicates based on key, and remove the key column
              waterfallFormat$mutation <- factor(waterfallFormat$mutation, levels=hierarchy$mutation)
              waterfallFormat <- waterfallFormat[order(waterfallFormat$mutation),]
              waterfallFormat <- waterfallFormat[!duplicated(waterfallFormat$key),]
              waterfallFormat[,key:=NULL]
              
              # print status message
              if(verbose){
                  memo <- paste("Removed", rowCountOrig - nrow(waterfallFormat),
                                "rows from the data which harbored duplicate",
                                "genomic locations")
                  message(memo)
              }
              
              # convert appropriate columns to factor
              waterfallFormat$sample <- factor(waterfallFormat$sample)
              
              return(waterfallFormat)
          })

#' @rdname toMutSpectra-methods
#' @aliases toMutSpectra,VEP
#' @importFrom BSgenome getSeq
#' @noRd
setMethod(f="toMutSpectra",
          signature="VEP",
          definition=function(object, BSgenome, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Converting", class(object),
                                "to expected MutSpectra format")
                  message(memo)
              }
              
              # grab the BSgenome
              if(is.null(BSgenome)){
                  if(verboase){
                      memo <- paste("Looking for correct genome for reference base annotation.")
                      messsage(memo)
                  }
                  
                  # look for assembly version in header
                  header <- object@vepObject@header
                  header <- header$Info[grepl("assembly", header$Info)]
                  header <- regmatches(header,regexpr("\\w(\\d)+", header))
                  if(length(header) != 1){
                      memo <- paste("Unable to infer assembly from VEP header")
                      stop(memo) 
                  }
                  
                  # determine if a genome is available
                  availableGenomes <- BSgenome::available.genomes()
                  availableGenomes <- availableGenomes[grepl(header, availableGenomes)]
                  if(length(availableGenomes) == 0){
                      memo <- paste("Could not find a compatible BSgenome for", toString(header),
                                    "Please specify the bioconductor BSgenome to annotate references bases!")
                      stop(memo)
                  }
                  
                  # determine if the available genome in an installed package
                  installedGenomes <- BSgenome::installed.genomes()
                  installedGenomes <- installedGenomes[installedGenomes == availableGenomes]
                  if(legnth(installedGenomes) == 0){
                      memo <- paste("The BSgenome", toString(availableGenomes), "is available",
                                    "but is not installed! Please install", toString(availableGenomes),
                                    "via bioconductor!")
                      stop(memo)
                  }
              }
              
              # grab the sample, mutation, position columns
              sample <- object@vepObject@sample
              variantAllele <- object@vepObject@mutation[,"Allele"]
              position <- object@vepObject@position[,"Location"]
              positionSplit <- lapply(as.character(position$Location), strsplit, ":", fixed=TRUE)
              chr <- unlist(lapply(positionSplit, function(x) x[[1]][1]))
              coord <- unlist(lapply(positionSplit, function(x) x[[1]][2]))
              coord <- lapply(as.character(tmp), strsplit, "-", fixed=TRUE)
              
              # combine all columns into a consistent format
              mutSpectraFormat <- cbind(sample, position, refAllele, variantAllele)
              colnames(mutSpectraFormat) <- c("sample", "position", "refAllele", "variantAllele")
              
              # unique, to make sure no duplicate variants exist to throw off the counts
              rowCountOrig <- nrow(mutSpectraFormat)
              mutSpectraFormat <- unique(mutSpectraFormat)
              
              # print status message
              if(verbose){
                  memo <- paste("Removed", rowCountOrig - nrow(mutSpectraFormat),
                                "rows from the data which harbored duplicate",
                                "genomic locations")
                  message(memo)
                  memo <- paste("Annotating reference bases")
                  message(memo)
              }
              tmp@vepObject@position
              

              # convert appropriate columns to factor
              mutSpectraFormat$sample <- factor(mutSpectraFormat$sample)
              
              return(mutSpectraFormat)
          })