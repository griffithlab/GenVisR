################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Class VarScanFormat
#' 
#' An S4 class acting as a container for VarScanFormat.
#' @name VarScanFormat-class
#' @rdname VarScanFormat-class
#' @slot path Character string specifying the path of the VarScan file read in.
#' @exportClass VarScanFormat
#' @include VarScanFormat_Virtual-class.R
#' @import methods
setClass("VarScanFormat",
         representation=representation(path="character"), 
         contains="VarScanFormat_Virtual",
         validity = function(object) {
             head(object@varscan)
             ## Expected varscan column names
             cnames <- c("chrom", "position", "ref", "var",
                         "normal_reads1", "normal_reads2", "normal_var_freq",
                         "normal_gt", "tumor_reads1", "tumor_reads2", "tumor_var_freq",
                         "tumor_gt", "somatic_status", "variant_p_value",
                         "somatic_p_value", "tumor_reads1_plus", "tumor_reads1_minus",
                         "tumor_reads2_plus", "tumor_reads2_minus",
                         "normal_reads1_plus", "normal_reads1_minus",
                         "normal_reads2_plus", "normal_reads2_minus", "sample")
             
             ## Check to see if there is any data after the filtering steps for varscan
             if (nrow(object@varscan) == 0) {
                 stop("No varscan data can be found after filtering based on normal VAF and Germline/LOH somatic_status")
             }
             
             ## Check the column names to see if there is the appropriate input
             varscan_column_names <- colnames(object@varscan)
             num <- which(!varscan_column_names%in%cnames)
             if (length(num) > 0 & length(varscan_column_names) == length(cnames)) {
                 stop("Column names of varscan input are not what is expected. Please
                      refer to 
                      http://varscan.sourceforge.net/somatic-calling.html#somatic-output 
                      for appropriate column names.")
             }
             if (length(num) > 0 & length(varscan_column_names) != length(cnames)) {
                 stop("Number of columns in varscan input are not what is expected. 23
                      columns are expected. Please refer to 
                      http://varscan.sourceforge.net/somatic-calling.html#somatic-output
                      for appropriate columns and column names.")
             }
             
             ## Check to see if the VAF columns are proportion as opposed to percentage
             ## Function requires input in percentages and will convert percentage to proportion
             tumor_per <- any(grepl("%", object@varscan$tumor_var_freq) == TRUE)
             normal_per <- any(grepl("%", object@varscan$normal_var_freq) == TRUE)
             if (tumor_per == TRUE | normal_per == TRUE) {
                 stop("Make sure the tumor/normal VAF column is in percentage and not proportion. 
                      (i.e. 75.00% as opposed to 0.75).")
             }
             
             ## Check to see if the VAF provided are somatic or not
             if (any(object@varscan$tumor_var_freq>1 | object@varscan$normal_var_freq >1)) {
                 stop("Detected values in either the normal or tumor variant ",
                      "allele fraction columns above 1. Values supplied should ",
                      "be a proportion between 0-1!")
             }
             if (any(object@varscan$normal_var_freq<0.4 | object@varscan$normal_var_freq>0.6)) {
                 stop("Detected values with a variant allele fraction either ",
                         "above .6 or below .4 in the normal. Please ensure ",
                         "variants supplied are heterozygous in the normal!. 
                         Make sure to remove coordinates with normal VAF > 0.6 or < 0.4")
             }
             
             return(TRUE)
         }
)

#' Constructor for the VarScanFormat container class.
#' 
#' @name VarScanFormat
#' @rdname VarScanFormat-class
#' @param path String specifying the path to a VarScan file.
#' @param verbose Boolean specifying if progress should be reported while reading
#' in the VarScan. file.
#' @seealso \code{\link{lohSpec}}
#' @importFrom data.table fread
#' @export
VarScanFormat <- function(path, verbose=FALSE) {
    ## Read in VarScan data
    varscanData <- suppressWarnings(fread(input=path, stringsAsFactors=FALSE,
                                                      verbose=verbose))
    
    ## Get the sample names
    sample <- varscanData[,which(colnames(varscanData)=="sample"), with=FALSE]

    ## Convert VAF percentages to VAF proportions
    varscanData$normal_var_freq <- round(as.numeric(as.character(gsub(pattern = "%", 
                                                                           replacement = "", varscanData$normal_var_freq)))/100, digits = 3)
    varscanData$tumor_var_freq <- round(as.numeric(as.character(gsub(pattern = "%", 
                                                                      replacement = "", varscanData$tumor_var_freq)))/100, digits = 3)
    
    ## Obtain coordinates that were called as germline or LOH by varscan
    varscanData <- varscanData[somatic_status == "Germline" | somatic_status == "LOH"]
    
    ## Remove coordinates with normal VAF > 0.6 or < 0.4
    varscanData <- varscanData[normal_var_freq<=0.6 &
                                         normal_var_freq>=0.4]
    
    ## Create the varscan object
    varscanObject <- new(Class="VarScanFormat", path=path, varscan=varscanData, sample=sample)
    return(varscanObject)
    
}

################################################################################
####################### Method function definitions ############################

#' @rdname getLohData-methods
#' @aliases getLohData
#' @noRd
#' @importFrom data.table data.table
setMethod(f="getLohData",
          signature="VarScanFormat",
          definition=function(object, verbose, ...) {

              ## Print status message
              if (verbose) {
                  message("Generating LOH dataset.")
              }
              
              ## Get the necessary columns from varscan output
              primaryData <- object@varscan[,c("chrom", "position", "tumor_var_freq", 
                                        "normal_var_freq", "sample"), 
                                     with=FALSE]
              colnames(primaryData) <- c("chromosome", "position", "tumor_var_freq", 
                                         "normal_var_freq", "sample")

              return(primaryData)
              
          })




