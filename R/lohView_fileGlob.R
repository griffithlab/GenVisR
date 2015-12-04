#' Grab data for lohView
#'
#' Look in the specified file path and grab data with the proper extension for
#' lohView
#' @name lohView_fileGlob
#' @param path character string specifying which directory contains 
#' the sample information stored as datasets with columns "chromosome", 
#' "position", "n_freq", "t_freq", and "sample" (required if x is not specified)
#' @param fileExt character string specifying the file extensions of files
#' @param step integer with the length of divisions (bp) in chromosomes
#' @return object of class data frame from data specified in path for lohView

lohView_fileGlob <- function(path, fileExt, step, gender)
{
    fileNames <- Sys.glob(paste0(path, '*.', fileExt))
    columnNames <- c("chromosome", "position", "n_freq", "t_freq", "sample")
    
    for (i in 1:length(fileNames))
    {
        data <- read.delim(fileNames[i])
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
    
    if (gender == TRUE) {
        dataset <- dataset[dataset$chromosome !="Y",]
    }
    if(gender == FALSE) {
        dataset <- dataset[dataset$chromosome != "X" & 
                               dataset$chromosome != "Y",]
    }
    
    return(dataset)
}