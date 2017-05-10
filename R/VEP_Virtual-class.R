#' Class VEP_Virtual
#' 
#' An S4 class to act as a virtual class for VEP version sub-classes.
#' @name VEP_Virtual-class
#' @rdname VEP_Virtual-class
#' @slot header data.table object holding header information.
#' @slot description data.table object holding column descriptions
#' @slot position data.table object holding genomic positions.
#' @slot mutation data.table object holding mutation status data.
#' @slot sample data.table object holding sample data.
#' @slot meta data.table object holding all other meta data.
#' @importClassesFrom data.table data.table
#' @import methods
setClass(
    Class="VEP_Virtual",
    representation=representation(header="data.table",
                                  description="data.table",
                                  position="data.table",
                                  mutation="data.table",
                                  sample="data.table",
                                  meta="data.table", "VIRTUAL")
)

#' @rdname parseDescription-methods
#' @aliases parseDescription,VEP_Virtual
#' @param object Object of class VEP_Virtual, needed for method dispatch on class.
#' @param vepHeader List of character vectors obtained from the vep header
#' @importFrom data.table setDT
setMethod(f="parseDescription",
          signature="VEP_Virtual",
          definition=function(object, vepHeader, ...){
              # anonymous function to grab only the column descriptions
              a <- function(x){
                  descriptionFieldIndex <- which(grepl("Extra column keys", x)) + 1
                  descriptionFieldIndex <- descriptionFieldIndex:length(x)
                  x <- x[descriptionFieldIndex]
                  return(x)
              }
              # obtain the column descriptions and clean up
              description<- lapply(vepHeader, a)
              description <- unique(unlist(description))
              description <- gsub("## ", "", description)
              
              # convert these results to a data.table after splitting into two columns
              description <- as.data.frame(do.call(rbind, strsplit(description, "\\s*:\\s*")))
              names(description) <- c("Name", "Description")
              data.table::setDT(description)
              
              # return the results
              return(description)
          })

#' @rdname parseHeader-methods
#' @aliases parseHeader,VEP_Virtual
#' @param object Object of class VEP_Virtual, needed for method dispatch on class.
#' @param vepHeader List of character vectors obtained from the vep header
setMethod(f="parseHeader",
          signature="VEP_Virtual",
          definition=function(object, vepHeader, ...){
              
              # anonymous function to grab only the column headers
              a <- function(x){
                  headerFieldIndex <- which(grepl("Extra column keys", x)) - 1
                  headerFieldIndex <- 1:headerFieldIndex
                  x <- x[headerFieldIndex]
                  return(x)
              }
              
              # obtain the column headers and clean up
              header <- lapply(vepHeader, a)
              header <- lapply(header, function(x) x[-which(grepl("Output produced at", x))])
              header <- unique(unlist(header))
              header <- gsub("## ", "", header)
              
              # convert these results to a data.table 
              header <- data.table::as.data.table(header)
              names(header) <- c("Info")
              
              # return the results
              return(header)
          })

#' @rdname parseExtra-methods
#' @aliases parseExtra,VEP_Virtual
#' @param object Object of class VEP_Virtual, needed for method dispatch on class.
#' @param vepData Object of class data.table holding vep annotated data
#' @import data.table
setMethod(f="parseExtra",
          signature="VEP_Virtual",
          definition=function(object, vepData, ...){
              
              # check that "extra" column is present if not return data as it was
              if(!any(colnames(vepData) %in% c("Extra"))){
                  return(vepData)
              }
              
              # Split fields in the "Extra" column of a VEP file into actual columns
              extraCol <- lapply(strsplit(as.character(vepData$Extra), ';', fixed=TRUE),
                                 function(x){
                                     res <- data.table::tstrsplit(x, '=', fixed=TRUE)
                                     setNames(as.list(res[[2]]), res[[1]])
                                     })
              extraCol <- rbindlist(extraCol, fill = T)
              
              # meged the Extra column back in
              vepData <- vepData[,Extra:=NULL]
              vepData <- cbind(vepData, extraCol)
              
              return(vepData)
          })

#' @rdname writeData-methods
#' @aliases writeData,VEP_Virtual
#' @param object Object of class VEP_Virtual.
#' @param file Character string specifying a file to send output to.
#' @param sep Delimiter used when writing output, defaults to "\t".
#' @importFrom data.table fwrite
setMethod(f="writeData",
          signature="VEP_Virtual",
          definition=function(object, file, sep, ...){
              data.table::fwrite(cbind(object@position, object@mutation, object@sample,
                                       object@meta), file=file, sep=sep, na=NA)
          })