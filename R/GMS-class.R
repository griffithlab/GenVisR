#' Class GMS
#' 
#' An S4 class for Genome Modeling System annotation files.
#' @name GMS-class
#' @rdname GMS-class
#' @slot path Character string specifying the paths of the GMS files read in.
#' @slot version Numeric value specifying the version of the GMS annotation files.
#' @slot gmsObject gms object which inherits from gms_Virtual class.
#' @exportClass GMS
#' @include GMS_Virtual-class.R
#' @import methods
setClass("GMS",
         representation=representation(path="character",
                                       version="numeric",
                                       gmsObject="GMS_Virtual"))

#' Initalizer method for the GMS container class
#' 
#' @name GMS
#' @rdname GMS-class
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
#' @importFrom data.table is.data.table
#' @importFrom data.table as.data.table
setMethod(
    f="initialize",
    signature="GMS",
    definition=function(.Object, path, data, version, verbose){
        cat("!!!!! GMS~Initalizer !!!!!\n")
        
        # assign the annotator version
        .Object@version <- version
        
        # if data is null read from path variable
        if(is.null(data)){
            # Grab all files and assign to slot
            gmsFiles <- Sys.glob(path)
            .Object@path <- gmsFiles
            
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
                # read data
                gmsData <- suppressWarnings(data.table::fread(input=x,
                                                              stringsAsFactors=TRUE,
                                                              verbose=verbose))
                # set sample if it's not already in the data table
                if(any(colnames(gmsData) %in% "sample")){
                    return(gmsData)
                } else {
                    gmsData$sample <- sampleName
                    return(gmsData)
                }
            }
            
            # aggregate data into a single data table if necessary
            if(length(gmsFiles) == 0){
                memo <- paste("No files found using:", path)
                stop(memo)
            } else {
                gmsData <- lapply(gmsFiles, a, verbose)
                gmsData <- data.table::rbindlist(gmsData)
            }
        } else if(is.data.table(data)){
            .Object@path <- "none"
            gmsData <- data
        } else {
            memo <- paste("data is not of class data.table,",
                          "attempting to coerce")
            warning(memo)
            .Object@path <- "none"
            gmsData <- data.table::as.data.table(data)
        }
        
        # assign the gmsData to it's slot
        if(version == 4){
            .Object@gmsObject <- GMS_v4(gmsData=gmsData)
        } else {
            memo <- paste("Currently only GMS version 4 is supported, make a",
                          "feature request on",
                          "https://github.com/griffithlab/GenVisR!")
            stop(memo)
        }
        
        return(.Object)
    })

#' Constructor for the GMS container class.
#' 
#' @name GMS
#' @rdname GMS-class
#' @param path String specifying the path to a GMS annotation file. Can accept
#' wildcards if multiple GMS annotation files exist (see details).
#' @param data data.table object storing a GMS annotation file. Overrides "path"
#' if specified.
#' @param version String specifying the version of the GMS files, Defaults to
#' version 4.
#' @param verbose Boolean specifying if progress should be reported while
#' reading in the GMS files.
#' @details When specifying a path to a GMS annotation file the option exist to
#' either specify the full path to an annotation file or to use wildcards to
#' specify multiple files. When specifying a full path the initalizer will check
#' if a column named "sample" containg the relevant sample for each row exists.
#' If such a column is not found the initalizer will assume this file
#' corresponds to only one sample and populate a sample column accordingly.
#' Alternatively if multiple files are specified at once using a wildcard, the
#' initalizer will aggregate all the files and use the file names minus any
#' extension top populate sample names.
#' The version defaults to 4 which is the default value of the GMS annotator.
#' This value will need to be changed only if files were created using a
#' different GMS annotator version.
#' @export
GMS <- function(path, data=NULL, version=4, verbose=FALSE){
    cat("!!!!! GMS~Constructor !!!!!\n")
    new("GMS", path=path, data=data, version=version, verbose=verbose)}

#' @rdname getPosition-methods
#' @aliases getPosition,GMS
setMethod(f="getPosition",
          signature="GMS",
          definition=function(object, ...){
              positions <- object@gmsObject@position
              return(positions)
          })

#' @rdname getMutation-methods
#' @aliases getMutation,GMS
setMethod(f="getMutation",
          signature="GMS",
          definition=function(object, ...){
              mutations <- object@gmsObject@mutation
              return(mutations)
          })

#' @rdname getSample-methods
#' @aliases getSample,GMS
setMethod(f="getSample",
          signature="GMS",
          definition=function(object, ...){
              sample <- object@gmsObject@sample
              return(sample)
          })

#' @rdname getMeta-methods
#' @aliases getMeta,GMS
setMethod(
    f="getMeta",
    signature="GMS",
    definition=function(object, ...){
        return(object@gmsObject@meta)
    }
)

#' @rdname toWaterfall-methods
#' @aliases toWaterfall,GMS
#' @noRd
setMethod(f="toWaterfall",
          signature="GMS",
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
              sample <- object@gmsObject@sample
              mutation <- object@gmsObject@mutation[,"trv_type"]
              gene <- object@gmsObject@meta[,"gene_name"]
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
                                    "in the meta slot of the gmsObject! Valid",
                                    "names are:", colnames(getMeta(object)))
                      warning(memo)
                      next
                  } else {
                      label <- object@gmsObject@meta[,labelColumn]
                  }
              }
              
              # combine all columns into a consistent format
              waterfallFormat <- cbind(sample, gene, mutation, label)
              colnames(waterfallFormat) <- c("sample", "gene", "mutation", "label")
              
              # make a temporary ID column for genomic features to collapse on
              # this will ensure the mutation burden/frequency plot will be accurate
              waterfallFormat$key <- paste0(object@gmsObject@position$chromosome_name, ":",
                                            object@gmsObject@position$start, ":",
                                            object@gmsObject@position$stop, ":",
                                            object@gmsObject@mutation$reference, ":",
                                            object@gmsObject@mutation$variant, ":",
                                            object@gmsObject@sample$sample)
              rowCountOrig <- nrow(waterfallFormat)
              
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

#' @rdname setMutationHierarchy-methods
#' @aliases setMutationHierarchy,GMS
#' @noRd
#' @importFrom data.table data.table
#' @importFrom data.table setDT
#' @importFrom data.table rbindlist
setMethod(f="setMutationHierarchy",
          signature="GMS",
          definition=function(object, mutationHierarchy, verbose, ...){
              # set the mutation hierarchy if a custom hierarchy is unspecified
              if(is.null(mutationHierarchy)){
                  mutationHierarchy$mutation <- c("nonsense", "frame_shift_del",
                                                  "frame_shift_ins", "splice_site_del",
                                                  "splice_site_ins", "splice_site",
                                                  "nonstop", "in_frame_del", "in_frame_ins",
                                                  "missense", "splice_region_del",
                                                  "splice_region_ins", "splice_region",
                                                  "5_prime_flanking_region", 
                                                  "3_prime_flanking_region", 
                                                  "3_prime_untranslated_region",
                                                  "5_prime_unstranslated_region",
                                                  "rna", "intronic", "silent")
                  
                  mutationHierarchy$color <- c('#4f00A8', '#A80100', '#CF5A59',
                                               '#A80079', '#BC2D94', '#CF59AE',
                                               '#000000', '#006666', '#00A8A8',
                                               '#009933', '#ace7b9', '#cdf0d5',
                                               '#59CF74', '#002AA8', '#5977CF',
                                               '#F37812', '#F2B079', '#888811',
                                               '#FDF31C', '#8C8C8C')
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
              if(!all(object@gmsObject@mutation$trv_type %in% mutationHierarchy$mutation)){
                  missingMutations <- unique(object@gmsObject@mutation$trv_type[!object@gmsObject@mutation$trv_type %in% mutationHierarchy$mutation])
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