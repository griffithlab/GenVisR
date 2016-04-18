#' Grab data for lohSpec
#'
#' Look in the specified file path and grab data with the proper extension for
#' lohSpec
#' @name lohSpec_fileGlob
#' @param path character string specifying which directory contains
#' the sample information stored as datasets with columns "chromosome",
#' "position", "n_vaf", "t_vaf", and "sample" (required if x is not specified)
#' @param fileExt character string specifying the file extensions of files
#' @param step integer with the length of divisions (bp) in chromosomes
#' @param gender vector of length equal to the number of samples, consisting of
#' elements from the set {"M", "F"}
#' @param window_size Integer value specifying the size of the window in base
#' pairs in which to calculate the mean Loss of Heterozygosity.
#' @return object of class data frame from data specified in path for lohSpec
#' @importFrom utils read.delim

lohSpec_fileGlob <- function(path, fileExt, step, window_size, gender)
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
        data <- utils::read.delim(fileNames[i])
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

        if (!exists("dataset", inherits=FALSE))
        {
            dataset <- data
        } else if(exists("dataset", inherits=FALSE)) {
            temp <- data
            dataset <- rbind(dataset, temp)
            rm(temp)
        }

        rm(data)
    }
    ## Solves problem of user removing loh values for any autosome
    all_lev <- unique(droplevels(dataset$chromosome))
    sample <- unique(dataset$sample)
    total <- data.frame()
        for (r in 1:length(sample)) {
            df <- dataset[dataset$sample==sample[r],]
            dflevels <- unique(droplevels(df$chromosome))
            chrDiff <- setdiff(all_lev, dflevels)
            if (length(chrDiff) >= 1) {
                for (i in 1:length(chrDiff)) {
                    if(is.null(gender)==TRUE) {
                        d1 <- c(as.numeric(chrDiff[i]), step, 50, 50, 
                                as.character(sample))
                        d2 <- c(as.numeric(chrDiff[i]), 
                                step + window_size, 50, 50, 
                                as.character(sample))
                        df <- data.frame(rbind(df, d1, d2))
                    }
                    if(is.null(gender)==FALSE) {
                        d1 <- c(as.character(chrDiff[i]), step, 50, 50, 
                                as.character(sample[r]), 
                                as.character(gender[r]))
                        
                        d2 <- c(as.character(chrDiff[i]),step + window_size, 
                                50, 50, 
                                as.character(sample[r]), 
                                as.character(gender[r]))
                        df <- data.frame(rbind(df, d1, d2))
                    }
                }
            }
            total <- rbind(total, df)
        }
    return(total)
}

