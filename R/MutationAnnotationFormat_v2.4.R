setClass("MutationAnnotationFormat_v2.4",
         contains="MutationAnnotationFormat",
         validity=function(object){
             cat("********** MutationAnnotationFormat_v2.4: inspector **********")
             expected_col <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
                               "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1",
                               "Tumor_Seq_Allele2", "Tumor_Sample_Barcode")
             if(!all(expected_col %in% colnames(object@data1))){
                 memo <- paste("Missing the following required columns:",
                               toString(expected_col[!expected_col %in% colnames(object@data1)]))
                 stop(memo)
             }else{}
             return(TRUE)
         }
)

setMethod(
    f="initialize",
    signature="MutationAnnotationFormat_v2.4",
    definition=function(.Object, path, version, data){
        .Object@data1 <- data
        .Object@path <- path
        .Object@version <- version
        return(.Object)
    }
)