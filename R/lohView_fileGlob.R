#' Grab data for lohView
#'
#' Look in the specified file path and grab data with the proper extension for
#' lohView
#' @name lohView_fileGlob
#' @param path character string specifying which directory contains 
#' the sample information stored as datasets with columns "chromosome", 
#' "position", "n_vaf", "t_vaf", and "sample" (required if x is not specified)
#' @param fileExt character string specifying the file extensions of files
#' @param step integer with the length of divisions (bp) in chromosomes
#' @return object of class data frame from data specified in path for lohView

lohView_fileGlob <- function(path, fileExt, step, gender)
{
    # Obtain file names with the tumor and normal vaf values
    fileNames <- Sys.glob(paste0(path, '*.', fileExt))
    # Determine the column names of the dataset
    if (is.null(gender) == FALSE) {
        columnNames <- c("chromosome", "position", "n_vaf", "t_vaf", "sample", 
                         "gender")
    }
    if(is.null(gender) == TRUE) {
        columnNames <- c("chromosome", "position", "n_vaf", "t_vaf", "sample")
    }
    
    # Extract raw t_vaf and n_vaf values and merge the dataset
    for (i in 1:length(fileNames))
    {
        data <- read.delim(fileNames[i])
        if (is.null(gender) == FALSE) {
            data <- data[data$chromosome !="Y",]
            if (is.null(data$gender) == TRUE) {
                data$gender <- gender[i]
            }
        }
        if(is.null(gender) == TRUE) {
            data <- data[data$chromosome != "X" & 
                                   data$chromosome != "Y",]
        }
        if(!all(columnNames %in% colnames(data)))
        {
            memo <- paste0("Did not detect all of the required columns in the",
                           " following file:", fileNames[i], "... skipping")
            warning(memo)
            next
        }

        if (!exists("dataset"))
        {
            dataset <- data
        } else if(exists("dataset")) {
            temp <- data
            dataset <- rbind(dataset, temp)
            rm(temp)
        }
        
        rm(data)
    }
    
    tail(dataset)
    return(dataset)
}