################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Class MutationAnnotationFormat
#' 
#' An S4 class acting as a container for MutationAnnotationFormat version sub-classes.
#' @name MutationAnnotationFormat-class
#' @rdname MutationAnnotationFormat-class
#' @slot path Character string specifying the path of the MAF file read in.
#' @slot version Numeric value specifying the version of the MAF file.
#' @slot mafObject MutationAnnotationFormat object which inherits from
#' MutationAnnotationFormat_Virtual class.
#' @exportClass MutationAnnotationFormat
#' @include MutationAnnotationFormat_Virtual-class.R
#' @import methods
setClass("MutationAnnotationFormat",
         representation=representation(path="character",
                                       version="numeric",
                                       mafObject="MutationAnnotationFormat_Virtual")
)

#' Constructor for the MutationAnnotationFormat container class.
#' 
#' @name MutationAnnotationFormat
#' @rdname MutationAnnotationFormat-class
#' @param path String specifying the path to a MAF file.
#' @param version String specifying the version of the MAF file, if set to auto
#' the version will be obtained from the header in the MAF file.
#' @param verbose Boolean specifying if progress should be reported while reading
#' in the MAF file.
#' @seealso \code{\link{Waterfall}}
#' @seealso \code{\link{MutSpectra}}
#' @importFrom data.table fread
#' @export
MutationAnnotationFormat <- function(path, version="auto", verbose=FALSE){
    mafData <- suppressWarnings(data.table::fread(input=path,
                                                  stringsAsFactors=TRUE,
                                                  verbose=verbose))
    
    # grab the maf version
    if(toupper(version) == "AUTO"){
        # read the version
        mafVersion <- readLines(con=path, n=50)
        mafVersion <- mafVersion[which(grepl("^#version", mafVersion))]
        mafVersion <- as.numeric(as.character(gsub("#version\\W", "", mafVersion)))
        if(length(mafVersion) == 0) stop("Unable to infer the maf Version from header, please specify")
    } else {
        mafVersion <- version
    }
    
    ##### Obtain the appropriate MutationAnnotationFormat sub-class ######
    if(mafVersion == 1.0){
        mafObject <- MutationAnnotationFormat_v1.0(mafData=mafData)
    }else if(mafVersion == 2.0){
        mafObject <- MutationAnnotationFormat_v2.0(mafData=mafData)
    }else if(mafVersion == 2.1){
        mafObject <- MutationAnnotationFormat_v2.1(mafData=mafData)
    }else if(mafVersion == 2.2){
        mafObject <- MutationAnnotationFormat_v2.2(mafData=mafData)
    }else if(mafVersion == 2.3){
        mafObject <- MutationAnnotationFormat_v2.3(mafData=mafData)
    }else if(mafVersion == 2.4 ){
        mafObject <- MutationAnnotationFormat_v2.4(mafData=mafData)
    }else{
        memo <- paste("The maf version:", toString(mafVersion),
                      "is currently unsupported. Make a feature request on",
                      "https://github.com/griffithlab/GenVisR!")
        stop(memo)
    }

    new("MutationAnnotationFormat", path=path, version=mafVersion, mafObject=mafObject)
}

################################################################################
###################### Accessor function definitions ###########################

#' @rdname getVersion-methods
#' @aliases getVersion
setMethod(f="getVersion",
          signature="MutationAnnotationFormat",
          definition=function(object, ...){
              version <- object@version
              return(version)
          })

#' @rdname getPath-methods
#' @aliases getPath
setMethod(f="getPath",
          signature="MutationAnnotationFormat",
          definition=function(object, ...){
              path <- object@path
              return(path)
          })

#' @rdname getPosition-methods
#' @aliases getPosition
setMethod(f="getPosition",
          signature="MutationAnnotationFormat",
          definition=function(object, ...){
              positions <- getPosition(object@mafObject)
              return(positions)
          })

#' @rdname getMutation-methods
#' @aliases getMutation
setMethod(f="getMutation",
          signature="MutationAnnotationFormat",
          definition=function(object, ...){
              mutations <- getMutation(object@mafObject)
              return(mutations)
          })

#' @rdname getSample-methods
#' @aliases getSample
setMethod(f="getSample",
          signature="MutationAnnotationFormat",
          definition=function(object, ...){
              sample <- getSample(object@mafObject)
              return(sample)
          })

#' @rdname getMeta-methods
#' @aliases getMeta
setMethod(f="getMeta",
          signature="MutationAnnotationFormat",
          definition=function(object, ...){
              meta <- getMeta(object@mafObject)
              return(meta)
          })

################################################################################
####################### Method function definitions ############################

#' @rdname toWaterfall-methods
#' @aliases toWaterfall
#' @noRd
setMethod(f="toWaterfall",
          signature="MutationAnnotationFormat",
          definition=function(object, hierarchy, labelColumn, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Converting", class(object),
                                "to expected waterfall format")
                  message(memo)
              }
              
              # grab the sample, mutation, gene columns and set a label
              sample <- getSample(object)
              mutation <- getMutation(object)[,"Variant_Classification"]
              gene <- getMeta(object)[,"Hugo_Symbol"]
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
              
              # make a temporary ID column for genomic features to collapse on
              # this will ensure the mutation burden/frequency plot will be accurate
              waterfallFormat$key <- paste0(getPosition(object)$Chromosome, ":",
                                            getPosition(object)$Start_Position, ":",
                                            getPosition(object)$End_Position, ":",
                                            getPosition(object)$Reference_Allele, ":",
                                            getPosition(object)$Tumor_Seq_Allele1, ":",
                                            getPosition(object)$Tumor_Seq_Allele2, ":",
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
#' @importFrom grDevices colors
setMethod(f="setMutationHierarchy",
          signature="MutationAnnotationFormat",
          definition=function(object, mutationHierarchy, verbose, ...){
              
              # set the mutation hierarchy if a custom hierarchy is unspecified
              if(is.null(mutationHierarchy)){
                  mutationHierarchy$mutation <- c("Nonsense_Mutation", "Frame_Shift_Ins",
                                                  "Frame_Shift_Del", "Translation_Start_Site",
                                                  "Splice_Site", "Nonstop_Mutation",
                                                  "In_Frame_Ins", "In_Frame_Del",
                                                  "Missense_Mutation", "5\'Flank",
                                                  "3\'Flank", "5\'UTR", "3\'UTR", "RNA",
                                                  "Intron", "IGR", "Silent",
                                                  "Targeted_Region")

                  mutationHierarchy$color <- c("grey", '#A80100', '#CF5A59',
                                               '#A80079', '#CF59AE', '#000000', 
                                               '#9159CF', '#4f00A8', '#59CF74',
                                               '#00A8A8', '#79F2F2', '#006666', 
                                               '#002AA8', '#5977CF', '#F37812',
                                               '#F2B079', '#888811', '#FDF31C')
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
              if(!all(object@mafObject@mutation$Variant_Classification %in% mutationHierarchy$mutation)){
                  missingMutations <- unique(object@mafObject@mutation$Variant_Classification[!object@mafObject@mutation$Variant_Classification %in% mutationHierarchy$mutation])
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
                  memo <- paste("The mutation type",toString(as.character(duplicateMut)),
                                "was duplicated in the supplied mutationHierarchy!")
                  warning(memo)
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
#' @param object Object of class MutationAnnotationFormat
#' @param verbose Boolean specifying if status messages should be reported
#' @importFrom data.table rbindlist
#' @noRd
setMethod(f="toMutSpectra",
          signature="MutationAnnotationFormat",
          definition=function(object, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Converting", class(object),
                                "to expected MutSpectra format")
                  message(memo)
              }
              
              # get an index of only the snvs
              snvIndex <- which(object@mafObject@meta$type == "SNP")
              
              # status message
              if(verbose){
                  memo <- paste("Removing", nrow(object@mafObject@sample)-length(snvIndex),
                                "entries which are not of class SNP!")
                  message(memo)
              }
              
              # grab the sample, mutation, position columns
              sample <- object@mafObject@sample[snvIndex]
              variantAllele1 <- object@mafObject@mutation[snvIndex,"Tumor_Seq_Allele1"]
              variantAllele2 <- object@mafObject@mutation[snvIndex,"Tumor_Seq_Allele2"]
              refAllele <- object@mafObject@mutation[snvIndex,"Reference_Allele"]
              chr <- object@mafObject@position[snvIndex,"Chromosome"]
              start <- object@mafObject@position[snvIndex,"Start_Position"]
              stop <- object@mafObject@position[snvIndex,"End_Position"]
              
              # combine all columns into a consistent format and remove duplicate variants
              mutSpectraFormat_Allele1 <- cbind(sample, chr, start, stop, refAllele, variantAllele1)
              mutSpectraFormat_Allele2 <- cbind(sample, chr, start, stop, refAllele, variantAllele2)
              colnames(mutSpectraFormat_Allele1) <- c("sample", "chromosome", "start", "stop", "refAllele", "variantAllele")
              colnames(mutSpectraFormat_Allele2) <- c("sample", "chromosome", "start", "stop", "refAllele", "variantAllele")
              mutSpectraFormat <- data.table::rbindlist(list(mutSpectraFormat_Allele1, mutSpectraFormat_Allele2))
              
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
#' @param object Object of class MutationAnnotationFormat
#' @param verbose Boolean specifying if status messages should be reported
#' @importFrom data.table rbindlist
#' @noRd
setMethod(f="toRainfall",
          signature="MutationAnnotationFormat",
          definition=function(object, verbose, ...){
              
              # print status message
              if(verbose){
                  memo <- paste("Converting", class(object),
                                "to expected Rainfall format")
                  message(memo)
              }
              
              # grab the sample, mutation, position columns
              sample <- object@mafObject@sample[snvIndex]
              variantAllele1 <- object@mafObject@mutation[snvIndex,"Tumor_Seq_Allele1"]
              variantAllele2 <- object@mafObject@mutation[snvIndex,"Tumor_Seq_Allele2"]
              refAllele <- object@mafObject@mutation[snvIndex,"Reference_Allele"]
              chr <- object@mafObject@position[snvIndex,"Chromosome"]
              start <- object@mafObject@position[snvIndex,"Start_Position"]
              stop <- object@mafObject@position[snvIndex,"End_Position"]
              
              # combine all columns into a consistent format and remove duplicate variants
              rainfallFormat_Allele1 <- cbind(sample, chr, start, stop, refAllele, variantAllele1)
              rainfallFormat_Allele2 <- cbind(sample, chr, start, stop, refAllele, variantAllele2)
              colnames(rainfallFormat_Allele1) <- c("sample", "chromosome", "start", "stop", "refAllele", "variantAllele")
              colnames(rainfallFormat_Allele2) <- c("sample", "chromosome", "start", "stop", "refAllele", "variantAllele")
              rainfallFormat <- data.table::rbindlist(list(rainfallFormat_Allele1, rainfallFormat_Allele2))
              
              # Remove cases where there is not change between reference and variant
              rainfallFormat$refAllele <- as.character(rainfallFormat$refAllele)
              rainfallFormat$variantAllele <- as.character(rainfallFormat$variantAllele)
              alleleMatchIndex <- rainfallFormat$refAllele == rainfallFormat$variantAllele
              rainfallFormat <- rainfallFormat[!alleleMatchIndex,]
              if(verbose){
                  memo <- paste("Removed", length(alleleMatchIndex), "entries",
                                "where the reference allele matched the tumor allele")
                  message(memo)
              }
              
              # unique, to make sure no duplicate variants exist to throw off the counts
              rowCountOrig <- nrow(rainfallFormat)
              rainfallFormat <- unique(rainfallFormat)
              
              # print status message
              if(verbose){
                  memo <- paste("Removed", rowCountOrig - nrow(rainfallFormat),
                                "rows from the data which harbored duplicate",
                                "genomic locations")
                  message(memo)
              }
              
              # make sure no duplicate genomic locations exist to throw off the mutation distance calculation
              dupCoordIndex <- duplicated(rainfallFormat[,c("sample", "chromosome", "start")])
              if(sum(dupCoordIndex) > 0){
                  rowCountOrig <- nrow(rainfallFormat)
                  rainfallFormat <- rainfallFormat[!dupCoordIndex,]
                  memo <- paste("Removed", rowCountOrig-nrow(rainfallFormat), "entries with identical",
                                "start coordinates")
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
              
              return(rainfallFormat)
          })