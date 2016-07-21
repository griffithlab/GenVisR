#' Class MutationAnnotationFormat
#' 
#' An S4 class to represent data in mutation annotation format
#' @name MutationAnnotationFormat-class
#' @rdname MutationAnnotationFormat-class
#' @slot path String specifying the path of the MAF file read in.
#' @slot version String specifying the version of the MAF file.
#' @slot position data.table object with columns "Chromosome", "Start_Position", "End_Position" and "Strand"
#' @slot mutation data.table object with columns "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"
#' @slot sample data.table object with columns "Tumor_Sample_Barcode"
#' @slot meta data.table object containing additional information
#' @exportClass MutationAnnotationFormat
#' @importFrom data.table data.table
#' @import methods
setClass("MutationAnnotationFormat",
         representation=representation(
                                       "path"="character",
                                       "version"="numeric",
                                       "position"="data.table",
                                       "mutation"="data.table",
                                       "sample"="data.table",
                                       "meta"="data.table"
         ),
         validity=function(object){
             cat("!!!!! MutationAnnotationFormat~Inspector !!!!!\n")
         #     expected_col <- c("Hugo_Symbol", "Chromosome", "Start_position", "End_position",
         #                       "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1",
         #                       "Tumor_Seq_Allele2", "Tumor_Sample_Barcode")
         #     if(!all(expected_col %in% colnames(object@data_file))){
         #         memo <- paste("Missing the following required columns:",
         #                       toString(expected_col[!expected_col %in% colnames(object@data_file)]))
         #         stop(memo)
         # }else{}
             return(TRUE)
         }
         
)

#' Initalizer method of MutationAnnotationFormat class
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

        ##### Create the appropriate child of parent MutationAnnotationFormat ######
        if(mafVersion == 1.0){
            mafObject <- new("MutationAnnotationFormat_v1.0", path=path,
                             version=mafVersion, data1=mafdata1)
        }else if(mafVersion == 2.0){
            mafObject <- new("MutationAnnotationFormat_v2.0", path=path,
                             version=mafVersion, data1=mafdata1)
        }else if(mafVersion == 2.1){
            mafObject <- new("MutationAnnotationFormat_v2.1", path=path,
                             version=mafVersion, data1=mafdata1)
        }else if(mafVersion == 2.3){
            mafObject <- new("MutationAnnotationFormat_v2.3", path=path,
                             version=mafVersion, data1=mafdata1)
        }else if(mafVersion == 2.4){
            mafObject <- MutationAnnotationFormat_v2.4(mafPath=path,
                                                       mafVersion=mafVersion,
                                                       mafData=mafData)
        }else{}
        browser()
        .Object@path <- getPath(mafObject)
        .Object@version <- getVersion(mafObject)
        .Object@data_file <- getDataFile(mafObject)
        .Object@postion <- getPosition(mafObject)
        
        validObject(.Object)
        return(.Object)
    }
)

#' Constructor for MutationAnnotationFormat
#' 
#' @name MutationAnnotationFormat
#' @rdname MutationAnnotationFormat-class
#' @param path String specifying the path to a MAF file.
#' @param version String specifying the version of the MAF file, if set to auto
#' the version will be grabbed from the header in the MAF file.
#' @param verbose Boolean specifying if progress should be reported while reading
#' in the MAF file.
#' @export
MutationAnnotationFormat <- function(path, version="auto", verbose=FALSE){
    cat("!!!!! MutationAnnotationFormat~Constructor !!!!!\n")
    new("MutationAnnotationFormat", path=path, version=version, verbose=verbose)
}
