setClass("MutationAnnotationFormat_v2.4",
         contains="MutationAnnotationFormat",
         validity=function(object){
             cat("!!!!! MutationAnnotationFormat_v2.4~Inspector !!!!!\n")
             expected_col <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
                               "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1",
                               "Tumor_Seq_Allele2", "Tumor_Sample_Barcode")
             if(!all(expected_col %in% colnames(object@data_file))){
                 memo <- paste("Missing the following required columns:",
                               toString(expected_col[!expected_col %in% colnames(object@data_file)]))
                 stop(memo)
             }else{}
             return(TRUE)
         }
)

setMethod(
    f="initialize",
    signature="MutationAnnotationFormat_v2.4",
    definition=function(.Object, path, version, data_file){
        cat("!!!!! MutationAnnotationVersion_v2.4~Initalizer !!!!!\n")
        .Object@data_file <- data_file
        .Object@path <- path
        .Object@version <- version
        validObject(.Object)
        return(.Object)
    }
)

MutationAnnotationFormat_v2.4 <- function(path, version, data_file){
    cat("!!!!! MutationAnnotationFormat_v2.4~Constructor !!!!!\n")
    new("MutationAnnotationFormat_v2.4", path=path, version=version, data_file=data_file)
}