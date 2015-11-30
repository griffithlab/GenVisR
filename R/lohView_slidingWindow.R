#' obtain LOH data
#' 
#' Obtain LOH heatmap on entire chromsomes from samples in a cohort
#' @name lohView_slidingWindow
#' @param path character string specifying which directory contains 
#' the sample information stored as datasets with columns "chromosome", 
#' "position", "n_freq", "t_freq", and "sample"
#' @param step integer with the length of divisions (bp) in chromosomes
#' @param window_size integer with the size of the sliding window (bp) to be 
#' applied
#' @param normal integer specifying the normal VAF frequency used in LOH 
#' calculations
#' @return object of class dataframe containing LOH data
#' @import dplyr
#' @examples
#' lohView_slidingWindow(path, step=500000, window_size=1000000, normal=50)
#' @export

lohView_slidingWindow <- function(path, step, window_size, normal) {
    setwd(path)
    fileNames <- Sys.glob("*.txt")
    for (i in 1:length(fileNames)) {
        data <- read.delim(fileNames[i])
        colnames(data) <- c("chromosome", "position", "n_freq", "tfreq", 
                            "sample")
        sample <- as.character(data$sample[i])
        chromosome <- as.character("Y")
        levels(data$chromosome) <- c(levels(data$chromosome), "Y")
        if (data$chromosome[nrow(data)] == "X") {
            new_data <- c(chromosome,step*.4,50,50,sample)
            new_data_2 <- c(chromosome,step*1.6,50,50,sample)
            new_data_3 <- c(chromosome,step*2.8,50,50,sample)
            new_data_4 <- c(chromosome,step*4,50,50,sample)
            data <- rbind(data, new_data, new_data_2, new_data_3, new_data_4)
        }
        if (!exists("dataset")) {
            dataset <- data
        }    
        else if(exists("dataset")) {
            temp <- data
            dataset <- rbind(dataset, temp)
            rm(temp)
        }
        rm(data)
    }
    loh_data <- dataset
    rm(dataset)
    colnames(loh_data) <- c("chromosome", "position", "n_freq", 
                            "t_freq", "sample")
    out <- split(loh_data, list(as.character(loh_data$chromosome),
                                      as.character(loh_data$sample)))

    window_data <- lohView_windowPosition(out, step, window_size)
    
    total <- data.frame()
    final_dataset <- list()
   for (i in 1:length(out)){
        final_dataset[[i]] <- plyr::adply(window_data[[i]], 1, lohView_lohCalc, 
                                          out[[i]], normal)
        
    }
   
    #loh_final_dataset <- lapply(window_data, function(x, out) plyr::adply(x, 1, lohView_lohCalc, out), out)
    
    
    loh_final_dataset <- plyr::ldply(final_dataset, data.frame)
    loh_dataset <- loh_final_dataset[,-c(4,5,6)]
    return(loh_dataset)
}

