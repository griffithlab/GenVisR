################################################################################
############################ Class Definition ##################################
setClass("MutationAnnotationFormat",
         representation=representation(
                                       "path"="character",
                                       "version"="numeric",
                                       "data_file"="data.table"
         ),
         validity=function(object){
             cat("!!!!! MutationAnnotationFormat~Inspector !!!!!\n")
             expected_col <- c("Hugo_Symbol", "Chromosome", "Start_position", "End_position",
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

################################################################################
############################ Initalizer ########################################
setMethod(
    f="initialize",
    signature="MutationAnnotationFormat",
    definition=function(.Object, path, version, verbose){
        cat("!!!!! MutationAnnotationFormat~Initalizer !!!!!\n")
        ##### Read in the file from the specified path ####
        # read the main data1
        mafdata1 <- suppressWarnings(data.table::fread(input=path,
                                                      stringsAsFactors=TRUE,
                                                      verbose=verbose))
        # grab the maf version
        if(version == "auto"){
            # read the version
            mafVersion <- readLines(con=path, n=50)
            mafVersion <- mafVersion[which(grepl("^#version", mafVersion))]
            mafVersion <- as.numeric(as.character(gsub("#version\\W", "", mafVersion)))
            if(length(mafVersion) == 0) stop("Unable to infer the maf Version from header, please specify via the parameter version")
        } else {
            mafVersion <- version
        }
        
        ##### Create the appropriate child of parent MutationAnnotationFormat ######
        if(mafVersion == 1.0){
            mafObject <- methods::new("MutationAnnotationFormat_v1.0", path=path,
                                      version=mafVersion, data1=mafdata1)
        }else if(mafVersion == 2.0){
            mafObject <- methods::new("MutationAnnotationFormat_v2.0", path=path,
                                      version=mafVersion, data1=mafdata1)
        }else if(mafVersion == 2.1){
            mafObject <- methods::new("MutationAnnotationFormat_v2.1", path=path,
                                      version=mafVersion, data1=mafdata1)
        }else if(mafVersion == 2.3){
            mafObject <- methods::new("MutationAnnotationFormat_v2.3", path=path,
                                      version=mafVersion, data1=mafdata1)
        }else if(mafVersion == 2.4){
            mafObject <- MutationAnnotationFormat_v2.4(path=path,
                                                       version=mafVersion,
                                                       data_file=mafdata1)
        }else{
            
        }
        
        mafObject <- as.MutationAnnotationFormat(mafObject)
        validObject(.Object)
        return(.Object)
    }
)

################################################################################
######################### Constructor ##########################################
MutationAnnotationFormat <- function(path, version="auto", verbose=FALSE){
    cat("!!!!! MutationAnnotationFormat~Constructor !!!!!\n")
    new("MutationAnnotationFormat", path=path, version=version, verbose=verbose)
}

################################################################################
################## additional Generics and Methods #############################

# setGeneric(
#     name="as.MutationAnnotationFormat",
#     def=function(object){standardGeneric("as.MutationAnnotationFormat")}
# )
# 
# setMethod(
#     f="as.MutationAnnotationFormat",
#     signature="MutationAnnotationFormat_v2.4",
#     definition=function(object){
#         colnames(object@data) <- c("Hugo_Symbol", "Chromosome", "Start_position", "End_position",
#                                     "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1",
#                                     "Tumor_Seq_Allele2", "Tumor_Sample_Barcode")
#         return(object)
#     }
# )

#MutationAnnotationFormat(path="~/Desktop/tcga_laml.maf")
