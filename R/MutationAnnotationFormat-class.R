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

#' Initalizer method for the MutationAnnotationFormat container class
#' 
#' @name MutationAnnotationFormat
#' @rdname MutationAnnotationFormat-class
#' @importFrom data.table fread
setMethod(
    f="initialize",
    signature="MutationAnnotationFormat",
    definition=function(.Object, path, version, verbose){
        
        cat("!!!!! MutationAnnotationFormat~Initalizer !!!!!\n")

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
        .Object@path <- path
        .Object@version <- mafVersion
        .Object@mafObject <- mafObject
        
        return(.Object)
    }
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
#' @export
MutationAnnotationFormat <- function(path, version="auto", verbose=FALSE){
    cat("!!!!! MutationAnnotationFormat~Constructor !!!!!\n")
    new("MutationAnnotationFormat", path=path, version=version, verbose=verbose)
}

#' @rdname getPosition-methods
#' @aliases getPosition,MutationAnnotationFormat
setMethod(f="getPosition",
          signature="MutationAnnotationFormat",
          definition=function(object, ...){
              positions <- object@mafObject@position
              return(positions)
          })

#' @rdname getMutation-methods
#' @aliases getMutation,MutationAnnotationFormat
setMethod(f="getMutation",
          signature="MutationAnnotationFormat",
          definition=function(object, ...){
              mutations <- object@mafObject@mutation
              return(mutations)
          })

#' @rdname getSample-methods
#' @aliases getSample,MutationAnnotationFormat
setMethod(f="getSample",
          signature="MutationAnnotationFormat",
          definition=function(object, ...){
              sample <- object@mafObject@sample
              return(sample)
          })

#' @rdname getMeta-methods
#' @aliases getMeta,MutationAnnotationFormat
setMethod(f="getMeta",
          signature="MutationAnnotationFormat",
          definition=function(object, ...){
              meta <- object@mafObject@meta
              return(meta)
          })

#' @rdname toWaterfall-methods
#' @aliases toWaterfall,MutationAnnotationFormat
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
              
              # grab the mutation hierarchy
              hierarchy <- hierarchy@MutationHierarchy
              
              # grab the sample, mutation, gene columns and set a label
              sample <- object@mafObject@sample
              mutation <- object@mafObject@mutation[,"Variant_Classification"]
              gene <- object@mafObject@meta[,"Hugo_Symbol"]
              label <- NA
              
              # if a label column exists and is proper overwrite the label variable
              if(!is.null(labelColumn)){
                  if(length(labelColumn) != 1) {
                      memo <- paste("Parameter \"labelColumn\" must be of length 1!",
                                    "Found length to be", length(labelColumn))
                      warning(memo)
                      next
                  } else if(labelColumn %in% colnames(object@mafObject@meta)){
                      memo <- paste("Did not find column:", labelColumn,
                                    "in the meta slot of the mafObject! Valid",
                                    "names are:", colnames(getMeta(object)))
                      warning(memo)
                      next
                  } else {
                      label <- object@mafObject@meta[,labelColumn]
                  }
              }
              
              # combine all columns into a consistent format
              waterfallFormat <- cbind(sample, gene, mutation, label)
              colnames(waterfallFormat) <- c("sample", "gene", "mutation", "label")
              
              # make a temporary ID column for genomic features to collapse on
              # this will ensure the mutation burden/frequency plot will be accurate
              waterfallFormat$key <- paste0(object@mafObject@position$Chromosome, ":",
                                            object@mafObject@position$Start_Position, ":",
                                            object@mafObject@position$End_Position, ":",
                                            object@mafObject@mutation$Reference_Allele, ":",
                                            object@mafObject@mutation$Tumor_Seq_Allele1, ":",
                                            object@mafObject@mutation$Tumor_Seq_Allele2, ":",
                                            object@mafObject@sample$sample)
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
#' @aliases setMutationHierarchy,MutationAnnotationFormat
#' @noRd
#' @importFrom data.table data.table
#' @importFrom data.table setDT
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

