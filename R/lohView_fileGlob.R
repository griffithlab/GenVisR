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

lohView_fileGlob <- function(path, fileExt, step)
{
    fileNames <- Sys.glob(paste0(path, '*', fileExt))
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
        
        sample <- as.character(data$sample[i])
        chromosome <- as.character("Y")
        levels(data$chromosome) <- c(levels(data$chromosome), "Y")
        
        if (data$chromosome[nrow(data)] == "X")
        {
            new_data <- c(chromosome,step*.4,50,50,sample)
            new_data_2 <- c(chromosome,step*1.6,50,50,sample)
            new_data_3 <- c(chromosome,step*2.8,50,50,sample)
            new_data_4 <- c(chromosome,step*4,50,50,sample)
            data <- rbind(data, new_data, new_data_2, new_data_3, new_data_4)
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
    
    return(dataset)
}