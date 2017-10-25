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
             head(object)
             ## Expected varscan column names
             cnames <- c("chrom", "position", "ref", "var",
                         "normal_reads1", "normal_reads2", "normal_var_freq",
                         "normal_gt", "tumor_reads1", "tumor_reads2", "tumor_var_freq",
                         "tumor_gt", "somatic_status", "variant_p_value",
                         "somatic_p_value", "tumor_reads1_plus", "tumor_reads1_minus",
                         "tumor_reads2_plus", "tumor_reads2_minus",
                         "normal_reads1_plus", "normal_reads1_minus",
                         "normal_reads2_plus", "normal_reads2_minus", "sample")
             
             ## Check the column names to see if there is the appropriate input
             varscan_column_names <- colnames(object@varscanData)
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
    browser()
    ## Read in VarScan data
    varscanData <- suppressWarnings(data.table::fread(input=path, 
                                                      stringsAsFactors=FALSE,
                                                      verbose=verbose))
    ## Put in sample name for now
    varscanData$sample <- "HCC1395"
    ## Get the sample names
    sample <- varscanData[,which(colnames(varscanData)=="sample"), with=FALSE]
    
    ## Create the varscan object
    varscanObject <- new("VarScanFormat", path=path, varscan=varscanData, sample=sample)
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
          definition=function(object, chr, verbose, ...) {
              ## Get the necessary columns from varscan output
              primaryData <- object@varscan[,c("chrom", "position", "tumor_var_freq", 
                                        "normal_var_freq", "sample"), 
                                     with=FALSE]
              
              ## Convert percentages to proportion
              primaryData$tumor_var_freq <- gsub("%", "", 
                                                 primaryData$tumor_var_freq)
              primaryData$normal_var_freq <- gsub("%", "", 
                                                  primaryData$normal_var_freq)
              primaryData$tumor_var_freq <- round(as.numeric(as.character(
                  primaryData$tumor_var_freq))/100, 
                                           digits = 3)
              primaryData$normal_var_freq <- round(as.numeric(as.character(
                  primaryData$normal_var_freq))/100, 
                                           digits = 3)
              
              ## Remove contigs, MT, and other unnecessary chromosomes
              if (is.null(chr)) {
                  chr <- c(as.character(seq(1:22)))
              }
              if (is.null(chr) == FALSE) {
                  primaryData <- primaryData[chrom %in% chr]
              }
              
              ## Print status message
              if (verbose) {
                  message("Generating LOH dataset.")
              }
              return(primaryData)
              
          })




