################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

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
#' @seealso \code{\link{Waterfall}}
#' @seealso \code{\link{MutSpectra}}
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
#' @importFrom data.table is.data.table
#' @importFrom data.table as.data.table
#' @export
GMS <- function(path, data=NULL, version=4, verbose=FALSE){
    
    # assign the annotator version
    version <- version
    
    # if data is null read from path variable
    if(is.null(data)){
        # Grab all files and assign to slot
        gmsFiles <- Sys.glob(path)
        path <- gmsFiles
        
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
        path <- "none"
        gmsData <- data
    } else {
        memo <- paste("data is not of class data.table,",
                      "attempting to coerce")
        warning(memo)
        path <- "none"
        gmsData <- data.table::as.data.table(data)
    }
    
    # assign the gmsData to it's slot
    if(version == 4){
        gmsObject <- GMS_v4(gmsData=gmsData)
    } else {
        memo <- paste("Currently only GMS version 4 is supported, make a",
                      "feature request on",
                      "https://github.com/griffithlab/GenVisR!")
        stop(memo)
    }
    
    new("GMS", path=path, gmsObject=gmsObject, version=version)
}

################################################################################
###################### Accessor function definitions ###########################

#' @rdname writeData-methods
#' @aliases writeData
setMethod(f="writeData",
          signature="GMS",
          definition=function(object, file, ...){
              writeData(object@gmsObject, file, sep="\t")
          })

#' @rdname getVersion-methods
#' @aliases getVersion
setMethod(f="getVersion",
          signature="GMS",
          definition=function(object, ...){
              version <- object@version
              return(version)
          })

#' @rdname getPath-methods
#' @aliases getPath
setMethod(f="getPath",
          signature="GMS",
          definition=function(object, ...){
              path <- object@path
              return(path)
          })

#' @rdname getPosition-methods
#' @aliases getPosition
setMethod(f="getPosition",
          signature="GMS",
          definition=function(object, ...){
              positions <- getPosition(object@gmsObject)
              return(positions)
          })

#' @rdname getMutation-methods
#' @aliases getMutation
setMethod(f="getMutation",
          signature="GMS",
          definition=function(object, ...){
              mutations <- getMutation(object@gmsObject)
              return(mutations)
          })

#' @rdname getSample-methods
#' @aliases getSample
setMethod(f="getSample",
          signature="GMS",
          definition=function(object, ...){
              sample <- getSample(object@gmsObject)
              return(sample)
          })

#' @rdname getMeta-methods
#' @aliases getMeta
setMethod(
    f="getMeta",
    signature="GMS",
    definition=function(object, ...){
        meta <- getMeta(object@gmsObject)
        return(meta)
    }
)

################################################################################
####################### Method function definitions ############################

#' @rdname toWaterfall-methods
#' @aliases toWaterfall
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
              
              # grab the sample, mutation, gene columns and set a label
              sample <- getSample(object)
              mutation <- getMutation(object)[,"trv_type"]
              gene <- getMeta(object)[,"gene_name"]
              label <- NA
              labelFlag <- TRUE
              
              # if a label column exists and is proper overwrite the label variable
              # if not change the flag
              if(!is.null(labelColumn)){
                  if(length(labelColumn) != 1) {
                      memo <- paste("Parameter \"labelColumn\" must be of length 1!",
                                    "Found length to be", length(labelColumn))
                      warning(memo)
                      labelFlag <- FALSE
                  }
                  
                  if(!labelColumn %in% colnames(getMeta(object))){
                      memo <- paste("Did not find column:", labelColumn,
                                    "in the meta slot of the vepObject! Valid",
                                    "names are:", toString(colnames(getMeta(object))))
                      warning(memo)
                      labelFlag <- FALSE
                  }
                  
                  if(labelFlag){
                      label <- getMeta(object)[,labelColumn, with=FALSE]
                  }
              }
              
              # combine all columns into a consistent format
              waterfallFormat <- cbind(sample, gene, mutation, label)
              colnames(waterfallFormat) <- c("sample", "gene", "mutation", "label")
              
              # remove any rows with NA values in gene or sample, this could mess up downstream functions
              rowCountOrig <- nrow(waterfallFormat)
              waterfallFormat <- waterfallFormat[!is.na(waterfallFormat$gene_name),]
              waterfallFormat <- waterfallFormat[!is.na(waterfallFormat$sample),]
              if(verbose & nrow(rowCountOrig) != nrow(waterfallFormat)){
                  memo <- paste("Removed", rowCountOrig-nrow(waterfallFormat), "rows which harbored NA values",
                                "in either the gene or sample column")
                  message(memo)
              }
              
              # make a temporary ID column for genomic features to collapse on
              # this will ensure the mutation burden/frequency plot will be accurate
              waterfallFormat$key <- paste0(getPosition(object)$chromosome_name, ":",
                                            getPosition(object)$start, ":",
                                            getPosition(object)$stop, ":",
                                            getPosition(object)$reference, ":",
                                            getPosition(object)$variant, ":",
                                            getSample(object)$sample)
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
#' @aliases setMutationHierarchy
#' @noRd
#' @importFrom data.table data.table
#' @importFrom data.table setDT
#' @importFrom data.table rbindlist
#' @importFrom grDevices colors
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
              if(!all(colnames(mutationHierarchy) %in% c("mutation", "color"))){
                  missingCol <- colnames(mutationHierarchy)[!c("mutation", "color") %in% colnames(mutationHierarchy)]
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
                  newCol <- grDevices::colors(distinct=TRUE)[!grepl("^gray", grDevices::colors(distinct=TRUE))]
                  tmp <- data.table::data.table("mutation"=missingMutations,
                                                "color"=sample(newCol, length(missingMutations)))
                  mutationHierarchy <- data.table::rbindlist(list(mutationHierarchy, tmp), use.names=TRUE, fill=TRUE)
              }
              
              # add in a pretty print mutation labels
              mutationHierarchy$label <- gsub("_", " ", mutationHierarchy$mutation)
              mutationHierarchy$label <-  gsub("'", "' ", mutationHierarchy$label)
              
              # check for duplicate mutations
              if(any(duplicated(mutationHierarchy$mutation))){
                  duplicateMut <- mutationHierarchy[duplicated(mutationHierarchy$mutation),"mutation"]
                  memo <- paste("The mutation type",toString(duplicateMut),
                                "was duplicated in the supplied mutationHierarchy!")
                  mutationHierarchy <- mutationHierarchy[!duplicated(mutationHierarchy$mutation),]
              }
              
              # ensure columns are of the proper type
              mutationHierarchy$color <- as.character(mutationHierarchy$color)
              mutationHierarchy$mutation <- as.character(mutationHierarchy$mutation)
              
              # print status message
              if(verbose){
                  memo <- paste("Setting the hierarchy of mutations from most",
                                "to least deleterious and mapping to colors:",
                                toString(mutationHierarchy$mutation))
                  message(memo)
              }
              
              return(mutationHierarchy)
          })

#' @rdname toMutSpectra-methods
#' @aliases toMutSpectra
#' @param object Object of class GMS
#' @param verbose Boolean specifying if status messages should be reported
#' @importFrom data.table rbindlist
#' @noRd
setMethod(f="toMutSpectra",
          signature="GMS",
          definition=function(object, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Converting", class(object),
                                "to expected MutSpectra format")
                  message(memo)
              }
              
              # get an index of only the snvs
              snvIndex <- which(object@gmsObject@meta$type == "SNP")
              
              # status message
              if(verbose){
                  memo <- paste("Removing", nrow(object@gmsObject@sample)-length(snvIndex),
                                "entries which are not of class SNP!")
                  message(memo)
              }
              
              # grab the sample, mutation, position columns
              sample <- object@gmsObject@sample[snvIndex]
              variantAllele <- object@gmsObject@mutation[snvIndex,"variant"]
              refAllele <- object@gmsObject@mutation[snvIndex,"reference"]
              chr <- object@gmsObject@position[snvIndex,"chromosome_name"]
              start <- object@gmsObject@position[snvIndex,"start"]
              stop <- object@gmsObject@position[snvIndex,"stop"]
              
              # combine all columns into a consistent format and remove duplicate variants
              mutSpectraFormat <- cbind(sample, chr, start, stop, refAllele, variantAllele)
              colnames(mutSpectraFormat) <- c("sample", "chromosome", "start", "stop", "refAllele", "variantAllele")
              
              # unique, to make sure no duplicate variants exist to throw off the counts
              rowCountOrig <- nrow(mutSpectraFormat)
              mutSpectraFormat <- unique(mutSpectraFormat)
              
              # print status message
              if(verbose){
                  memo <- paste("Removed", rowCountOrig - nrow(mutSpectraFormat),
                                "rows from the data which harbored duplicate",
                                "genomic locations")
                  message(memo)
              }
              
              # Remove cases where there is not change between reference and variant
              mutSpectraFormat$refAllele <- as.character(mutSpectraFormat$refAllele)
              mutSpectraFormat$variantAllele <- as.character(mutSpectraFormat$variantAllele)
              alleleMatchIndex <- mutSpectraFormat$refAllele == mutSpectraFormat$variantAllele
              mutSpectraFormat <- mutSpectraFormat[!alleleMatchIndex,]
              if(verbose){
                  memo <- paste("Removed", length(alleleMatchIndex), "entries",
                                "where the reference allele matched the tumor allele")
                  message(memo)
              }
              
              # convert appropriate columns to factor
              mutSpectraFormat$sample <- factor(mutSpectraFormat$sample)
              
              return(mutSpectraFormat)
          })

#' @rdname toRainfall-methods
#' @aliases toRainfall
#' @param object Object of class GMS
#' @param verbose Boolean specifying if status messages should be reported
#' @noRd
setMethod(f="toRainfall",
          signature="GMS",
          definition=function(object, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("converting", class(object), "to expected Rainfall format")
                  message(memo)
              }
              
              # grab the sample, position, and mutation columns
              sample <- getSample(object)
              chromosome <- getPosition(object)$chromosome_name
              start <- getPosition(object)$start
              stop <- getPosition(object)$stop
              refAllele <- as.character(getMutation(object)$reference)
              variantAllele <- as.character(getMutation(object)$variant)
              
              # combine all the relevant data into a single data table
              rainfallFormat <- cbind(sample, chromosome, start, stop, refAllele, variantAllele)
              
              # remove cases where a mutation does not exist
              rowCountOrig <- nrow(rainfallFormat)
              rainfallFormat <- rainfallFormat[rainfallFormat$refAllele != rainfallFormat$variantAllele,]
              if(rowCountOrig != nrow(rainfallFormat)){
                  memo <- paste("Removed", rowCountOrig - nrow(rainfallFormat),
                                "entries where the reference matches the variant")
                  warning(memo)
              }
              
              # remove mutations at duplicate genomic mutation as this could artifically increase
              # the density of mutations
              rowCountOrig <- nrow(rainfallFormat)
              rainfallFormat <- rainfallFormat[!duplicated(rainfallFormat[,c("sample", "chromosome", "start", "stop")]),]
              if(rowCountOrig != nrow(rainfallFormat)){
                  memo <- paste("Removed", rowCountOrig - nrow(rainfallFormat),
                                "entries with duplicate genomic positions")
                  warning(memo)
              }
              
              # create a flag column for where these entries came from
              rainfallFormat$origin <- 'mutation'
              
              # make sure everything is of the proper type
              rainfallFormat$sample <- factor(rainfallFormat$sample, levels=unique(rainfallFormat$sample))
              rainfallFormat$chromosome <- factor(rainfallFormat$chromosome, levels=unique(rainfallFormat$chromosome))
              rainfallFormat$start <- as.integer(rainfallFormat$start)
              rainfallFormat$stop <- as.integer(rainfallFormat$stop)
              rainfallFormat$refAllele <- factor(rainfallFormat$refAllele, levels=unique(rainfallFormat$refAllele))
              rainfallFormat$variantAllele <- factor(rainfallFormat$variantAllele, levels=unique(rainfallFormat$variantAllele))
              
              # return the standardized format
              return(rainfallFormat)
          })

#' @rdname toLolliplot-methods
#' @aliases toLolliplot
#' @param object Object of class GMS
#' @param verbose Boolean specifying if status messages should be reported
#' @noRd
setMethod(f="toLolliplot",
          signature="GMS",
          definition=function(object, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("converting", class(object), "to expected Lolliplot format")
                  message(memo)
              }
              
              # grab the sample, position, and mutation columns
              sample <- getSample(object)
              chromosome <- getPosition(object)$chromosome_name
              start <- getPosition(object)$start
              stop <- getPosition(object)$stop
              transcript <- getMeta(object)$transcript_name
              gene <- getMeta(object)$gene_name
              type <- getMutation(object)$trv_type
              aminoAcidChange <- getMeta(object)$amino_acid_change
              cPosition <- getMeta(object)$c_position
              
              # combine all the relevant data into a single data table
              lolliplotFormat <- cbind(sample, chromosome, start, stop, transcript, gene, type, aminoAcidChange, cPosition)
              
              return(lolliplotFormat)
              
              })
