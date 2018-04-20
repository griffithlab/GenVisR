#' Plot LOH data
#'
#' Construct a graphic visualizing Loss of Heterozygosity in a cohort
#' @name lohSpec
#' @param x object of class data frame with rows representing germline calls.
#' The data frame must contain columns with the following names "chromosome",
#' "position", "n_vaf", "t_vaf", "sample". required if path is set to NULL (see
#' details). vaf should range from 0-1.
#' @param y Object of class data frame with rows representing chromosome
#' boundaries for a genome assembly. The data frame must contain columns with
#' the following names "chromosome", "start", "end" (optional: see details).
#' @param genome Character string specifying a valid UCSC genome (see details).
#' @param gender Character vector of length equal to the number of samples,
#' consisting of elements from the set {"M", "F"}. Used to suppress the plotting
#' of allosomes where appropriate.
#' @param path Character string specifying the path to a directory containing
#' germline calls for each sample. Germline calls are expected to be stored as
#' tab-seperated files which contain the following column names "chromosome", 
#' "position", "n_vaf", "t_vaf", and "sample". required if x is set to null
#' (see details).
#' @param fileExt Character string specifying the file extensions of files
#' within the path specified. Required if argument is supplied to path
#' (see details).
#' @param step Integer value specifying the step size (i.e. the number of base
#' pairs to move the window). required when method is set to slide
#' (see details).
#' @param window_size Integer value specifying the size of the window in base
#' pairs in which to calculate the mean Loss of Heterozygosity (see details).
#' @param normal Numeric value within the range 0-1 specifying the expected
#' normal variant allele frequency to be used in Loss of Heterozygosity 
#' calculations. defaults to .50\%
#' @param colourScheme Character vector specifying the colour scale to use from
#' the viridis package. One of "viridis", "magma", "plasma", or "inferno".
#' @param plotLayer Valid ggpot2 layer to be added to the plot.
#' @param method character string specifying the approach to be used for 
#' displaying Loss of Heterozygosity, one of "tile" or "slide" (see details).
#' @param out Character vector specifying the the object to output, one of
#' "data", "grob", or "plot", defaults to "plot" (see returns).
#' @return One of the following, a list of dataframes containing data to be
#' plotted, a grob object, or a plot.
#' @details lohSpec is intended to plot the loss of heterozygosity (LOH) within
#' a sample. As such lohSpec expects input data to contain only LOH calls. Input
#' can be supplied as a single data frame given to the argument x with rows
#' containing germline calls and variables giving the chromosome, position, 
#' normal variant allele frequency, tumor variant allele frequency, and the
#' sample. In lieu of this format a series of .tsv files can be supplied via the 
#' path and fileExt arguments. If this method is choosen samples will be infered
#' from the file names. In both cases columns containing the variant allele
#' frequency for normal and tumor samples should range from 0-1.
#' Two methods exist to calculate and display LOH events. If the method is set
#' to "tile" mean LOH is calculated based on the window_size argument with
#' windows being placed next to each other. If the method is set to slide the
#' widnow will slide and calculate the LOH based on the step parameter.
#' In order to ensure the entire chromosome is plotted lohSpec requries the
#' location of chromosome boundaries for a given genome assembly. As a
#' convenience this information is available for the following genomes "hg19",
#' "hg38", "mm9", "mm10", "rn5" and can be tetrieved by supplying one of the
#' afore mentioned assemblies via the 'genome'paramter. If an argument is
#' supplied to the 'genome' parameter and is unrecognized a query to the UCSC
#' MySQL database will be attempted to obtain the required information. If
#' chromosome boundary locations are unavailable for a given assembly this
#' information can be supplied to the 'y' parameter which has priority over the
#' 'genome' parameter. 
#' @importFrom gtools mixedsort
#' @examples 
#' # plot loh within the example dataset
#' lohSpec(x=HCC1395_Germline)
#' @export

lohSpec <- function(x=NULL, path=NULL, fileExt=NULL, y=NULL, genome='hg19',
                    gender=NULL, step=1000000, window_size=2500000, 
                    normal=.50, colourScheme="inferno", plotLayer=NULL,
                    method="slide", out="plot")
{
    stop("This function has been deprecated in order to implement an object oriented programming style! Please use LohSpec() with a capital L instead!")
    
    # Grab data if necessary
    if(!is.null(path))
    {
        if(is.null(fileExt))
        {
            memo <- paste0("argument required to variable fileExt if argument ",
                           "is supplied to path")
            stop(memo)
        }
        x <- lohSpec_fileGlob(path=path, fileExt=fileExt, step=step, 
                              window_size=window_size, gender=gender)        
    }
    if (is.null(path)) {
        if (is.null(gender) == FALSE) {
            x <- x[x$chromosome !="Y",]
        }
        if(is.null(gender) == TRUE) {
            x <- x[(x$chromosome != "X" &
                        x$chromosome != "Y"),]
        }
    }
    
    # Data Quality Check
    input <- lohSpec_qual(x, y, genome)
    x <- input[[1]]
    y <- input[[2]]
    
    # Obtain dummy data for genome
    preloaded <- c('hg38', 'hg19', 'mm10', 'mm9', 'rn5')
    if(!is.null(y))
    {
        message("detected input to y, using supplied positions for chromosome
                boundaries")
        chr_pos <- y
    } else if(is.null(y) && any(genome == preloaded)) {
        message("genome specified is preloaded, retrieving data...")
        chr_pos <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == genome,]
        chr_pos <- multi_chrBound(chr_pos)
    } else {
        message("attempting to query UCSC sql database for chromosome
                positions")
        cyto_data <- suppressWarnings(multi_cytobandRet(genome))
        chr_pos <- multi_chrBound(cyto_data)
    }
    
    # Quality check for dummy data
    if(nrow(chr_pos) < 1)
    {
        memo <- paste0("Could not retrieve chromosome boundaries from",
                       " UCSC, please specify this information via ",
                       "the y paramter")
        stop(memo)
    }
    
    # Produce dataset with loh mean absolute differences 
    if (toupper(method) == 'SLIDE') {
        # Calculate loh via sliding window
        loh <- lohSpec_slidingWindow(loh_data=x, step, window_size, normal)
    }
    else if(toupper(method) == 'TILE') {
        # Calculate loh via tiled window
        ## Insert code
        loh <- lohSpec_tileWindow(loh_data=x, window_size, normal)
    }
    else {
        memo <- paste0("Did not recognize input to parameter method.", 
                       "Please specify one of \"Tile\" or \"Slide\".")
        stop(memo)
    }
    
    # set order of x axis labels in plot
    chromosomes <- gtools::mixedsort(as.character(unique(loh$chromosome)))
    
    # remove X and/or Y chromosomes
    if (is.null(gender) == FALSE) {
        chromosomes <- chromosomes[chromosomes != "Y"]
        chr_pos <- chr_pos[(chr_pos$chromosome != "chrY"),]
        loh <- loh[loh$chromosome != "Y",]
    }
    if (is.null(gender) == TRUE) {
        chromosomes <- chromosomes[chromosomes != "X" & chromosomes != "Y"]
        chr_pos <- chr_pos[(chr_pos$chromosome != "chrX" & 
                                chr_pos$chromosome != "chrY"),]
        loh <- loh[(loh$chromosome != "Y" & loh$chromosome != "X"),]
    }
    loh$chromosome <- factor(loh$chromosome, levels=chromosomes)
    chr_pos_levels <- gtools::mixedsort(as.character(unique(chr_pos$chromosome)))
    chr_pos$chromosome <- 
        factor(chr_pos$chromosome, levels=chr_pos_levels)
    
    # set order of y axis labels in plot
    samples <- gtools::mixedsort(as.character(unique(loh$sample)))
    loh$sample <- factor(loh$sample, levels=samples)
    
    #build  the plot
    loh_plot <- lohSpec_buildMain(loh, dummyData=chr_pos,
                                  colourScheme=colourScheme,
                                  plotLayer=plotLayer)
    
    # Decide what to output
    output <- multi_selectOut(data=loh, plot=loh_plot, draw=FALSE, out=out)
    return(output)
}

#' Plot LOH data
#'
#' Build a ggplot2 object displaying calculated LOH data
#' @name lohSpec_buildMain
#' @param x object of class dataframe with loh difference
#' calculations and column names "window_start", "window_stop", "chromosome",
#' "sample", and "loh_diff"
#' @param dummyData object of class dataframe with column names "chromosome",
#' "start", "end" specifying chromosome boundaries
#' @param colourScheme Character vector specifying the colour scale to use from
#' the viridis package. One of "viridis", "magma", "plasma", or "inferno".
#' @param plotLayer Valid ggplot2 layer to be added to the plot.
#' for the legend's parameters
#' @return object of class ggplot2
#' @import ggplot2
#' @importFrom viridis scale_fill_viridis

lohSpec_buildMain <- function(x, dummyData, colourScheme="inferno",
                              plotLayer=NULL)
{
    # define dummy data which will be chromosome boundaries, these are plotted
    # but are transparent and will not appear in the plot
    dummy_data <- geom_rect(data=dummyData, aes_string(xmin='start', xmax='end',
                                                       ymin=-1, ymax=1),alpha=0)
    # Define the main plot
    data <- geom_rect(data=x, aes_string(xmin='window_start',
                                         xmax='window_stop',
                                         ymin=-1,
                                         ymax=1, fill='loh_diff_avg'))
    
    
    # Define additional plot parameters
    facet <- facet_grid(sample ~ chromosome, scales="free", space="free")
    
    x_scale <- scale_x_continuous(expand = c(0, 0))
    y_scale <- scale_y_continuous(expand = c(0,0))
    
    lab_x <- xlab("Chromosome")
    lab_y <- ylab("Sample")
    
    # Define plot aesthetics
    BWscheme <- theme_bw()
    plotTheme <- theme(axis.ticks.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.text.y=element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank())
    
    # plot an additional layer if specified
    if(!is.null(plotLayer))
    {
        plotLayer <- plotLayer
    } else {
        plotLayer <- geom_blank()
    }
    
    # Define colour rame
    LOHgradient <- viridis::scale_fill_viridis("Avg. VAF Difference",
                                               option=colourScheme)
    
    # Build the plot
    tmp <- data.frame(x=0, y=0)
    p1 <- ggplot(data=tmp, aes(y=0)) + dummy_data + data + facet + x_scale + y_scale + 
        lab_x + lab_y + BWscheme + LOHgradient + plotTheme + plotLayer
    
    return(p1)
}

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
                            as.character(sample[r]))
                    d2 <- c(as.numeric(chrDiff[i]), 
                            step + window_size, 50, 50, 
                            as.character(sample[r]))
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

#' Calculate loh difference
#' 
#' Obtain LOH on an entire chromsomes from samples in a cohort
#' @name lohSpec_lohCalc
#' @param window_data object of class data frame with columns 
#' 'window_start' and 'window_stop
#' @param out object of class dataframe with columns 'chromosome', 
#' 'position', 'n_vaf', 't_vaf', and 'sample'
#' @param normal integer specifying the subtraction value from tumor VAF
#' @return object of class dataframe containing mean LOH difference calculations
#' and column names "window_start", "window_stop", "chromosome", "position", 
#' "n_vaf", "t_vaf", "sample", "loh_diff"

lohSpec_lohCalc <- function(window_data, out, normal)
    
{
    ## Set each list as a dataframe
    window <- data.frame(window_data)   
    sample_data <- data.frame(out)    
    sample_data$position <- as.numeric(as.character(sample_data$position))
    total <- data.frame()
    window_start <- vector()
    window_stop <- vector()
    num <- as.vector(1:nrow(window))
    
    df <- do.call("rbind", lapply(num, function(i){
        ## filter sample_data for loh calls whose positions fall in the
        ## window frame
        filtered_data <- 
            sample_data[sample_data$position >= window$window_start[i] &                           
                            sample_data$position <= window$window_stop[i],] 
        
        ## Perform the loh calculation to obtain avg loh in the window's frame
        loh_calc <- abs(as.numeric(as.character(filtered_data$t_vaf)) - normal)       
        loh_avg <- mean(loh_calc)
        window_start <- window$window_start[i]
        window_stop <- window$window_stop[i]
        if (is.na(loh_avg)==TRUE)    
        {            
            loh_avg <- NULL
            window_start <- NULL
            window_stop <- NULL
        }
        
        filtered_data$loh_diff_avg <- loh_avg
        filtered_data$window_start <- window_start
        filtered_data$window_stop <- window_stop
        return(filtered_data)
    }))
    
    return(df)
} 



#' Check input to lohSpec
#'
#' Perform data quality checks on input supplied to lohSpec
#' @name lohSpec_qual
#' @param x object of class data frame with columns 'chromosome', 'position',
#' 'n_vaf', 't_vaf', 'sample'
#' @param y object of class data frame with columns 'chromosome', 'start',
#' 'end' specifying chromosomal boundaries for a genome assembly
#' (required if genome is not specified)
#' @param genome character string specifying the genome assembly from which
#' input data is based
#' @return list of inputs passing basic quality controls

lohSpec_qual <- function(x, y, genome)
{
    # check input data to x
    if(!is.data.frame(x))
    {
        stop("Did not detect a data frame for input to x")
    }
    
    # check that correct columns are supplied in x
    x_col <- c('chromosome', 'position', 'n_vaf', 't_vaf', 'sample')
    if(!all(x_col %in% colnames(x)))
    {
        stop('Did not detect required column names in x, required columns are: '
             , paste0(x_col, sep="\t"))
    }
    
    # Check that values supplied in vaf columns are in the expected range
    if(any(x$nvaf > 1 | x$t_vaf > 1)) {
        memo <- paste("Detected values in either the normal or tumor variant",
                      "allele fraction columns above 1. Values supplied should",
                      "be a proportion between 0-1!")
        stop(memo)
    }
    
    if(any(x$n_vaf < .4 | x$n_vaf > .6)) {
        memo <- paste("Detected values with a variant allele fraction either",
                      "above .6 or below .4 in the normal. Please ensure",
                      "variants supplied are heterozygous in the normal!")
        warning(memo)
    }
    
    # Check chromosome column in x
    if(!all(grepl("^chr", x$chromosome)))
    {
        memo <- paste0("Did not detect the prefix chr in the chromosome column",
                       " of x... adding prefix")
        message(memo)
        x$chromosome <- paste0("chr", x$chromosome)
    } else if(all(grepl("^chr", x$chromosome))) {
        memo <- paste0("detected chr in the chromosome column of x...",
                       "proceeding")
        message(memo)
    } else {
        memo <- paste0("Detected unknown or mixed prefixes in the chromosome",
                       " column of x... should either be chr or none i.e. ",
                       "chr1 or 1")
        stop(memo)
    }
    
    # Check genome is acceptable name if y is not supplied
    if(is.null(y))
    {
        # Check that genome specified is not the ensembl name
        preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
        if(!any(genome == preloaded))
        {
            if(grepl("NCBI|GRC|RGSC|BROAD|BAYLOR|WUGSC",
                     genome, ignore.case=TRUE))
            {
                memo <- paste0("Detected a genome that does not appear to be,",
                               "in UCSC terms, please specify a genome in UCSC",
                               " terms to attempt query to UCSC mySQL databae.",
                               "Alternativly supply a value to y.")
                warning(memo)
            }
        }
    } else {
        
        # Check that y is a data frame
        if(!is.data.frame(y))
        {
            message("y is not a data frame, attempting to coerce")
            y <- as.data.frame(y)
        }
        
        # Check column names of y
        if(!all(c('chromosome', 'start', 'end') %in% colnames(y)))
        {
            memo <- paste0("Did not detect correct column names in y, missing",
                           "one of \"chromosome\", \"start\", \"end\"")
            stop(memo)
        }
        
        # Ensure that columns in data frame are of proper type
        y$chromosome <- as.character(y$chromosome)
        y$start <- as.integer(as.character(y$start))
        y$end <- as.integer(as.character(y$end))
    }
    
    return(list(x, y))
}

#' Obtain LOH data
#' 
#' Obtain LOH heatmap on entire chromsomes from samples in a cohort
#' @name lohSpec_slidingWindow
#' @param loh_data data frame with columns "chromosome", "position", "n_vaf",
#' "t_vaf", "sample" giving raw vaf calls for germline variants
#' @param step integer with the length of divisions (bp) in chromosomes
#' @param window_size integer with the size of the sliding window (bp) to be 
#' applied
#' @param normal integer specifying the normal VAF frequency used in LOH 
#' calculations
#' @return object of class dataframe containing LOH data
#' @importFrom plyr adply
#' @importFrom plyr ldply

lohSpec_slidingWindow <- function(loh_data, step, window_size, normal)
{     
    ## Obtain lists for each sample and chromosome
    out <- split(loh_data, list(as.character(loh_data$chromosome),
                                as.character(loh_data$sample)))
    
    ## Obtain the window position values 
    window_data <- lohSpec_windowPosition(out, step, window_size)
    
    total <- data.frame()
    final_dataset <- list()
    final_df <- list()
    loh_df <- list()
    ## Perform loh Calculations on each chromosome and sample within each window
    for (i in 1:length(out))
    {
        final_dataset[[i]] <- lohSpec_lohCalc(window_data[[i]], out[[i]], 
                                              normal)
        #final_dataset[[i]] <- plyr::ldply(final_df[[i]], data.frame)
    }
    
    ## Calculate avg loh for overlapping regions 
    df <- lohSpec_stepCalc(final_dataset, step = step) 
    
    ## Combine the lists into a single dataframe
    loh_dataset <- plyr::ldply(df, data.frame)
    colnames(loh_dataset) <- c("window_start", "window_stop", "chromosome", 
                               "sample", "loh_diff_avg")
    loh_dataset$loh_diff_avg <- loh_dataset$loh_diff_avg
    
    return(loh_dataset)
}

#' Obtain average loh within each step 
#' 
#' Calculate avverage LOH within each step
#' @name lohSpec_stepCalc
#' @param final_dataset object of class dataframe with columns 'window_start', 
#' 'window_stop', 'chromosome', 'position', 'n_vaf', 't_vaf', 'sample', and 
#' 'loh_diff_avg'
#' @param step integer with the length of divisions (bp) in chromosomes
#' @return list containing avg loh calculations for each step interval 

lohSpec_stepCalc <- function(final_dataset, step) {
    ## Remove irrelevant columns
    loh_df1 <- lapply(final_dataset, function(x) 
        x[!(names(x) %in% c("n_vaf", "t_vaf"))])
    
    step_final <- lapply(loh_df1, function(x){
        
        ## Obtain dataset with out repeating windowsloh
        loh_df <- x
        sample <- unique(loh_df$sample)
        chromosome <- unique(loh_df$chromosome)
        ## Obtain boundaries for each step sized window
        step_startMin <- min(as.integer(as.character(loh_df$window_start)))
        step_startMax <- max(as.integer(as.character(loh_df$window_start)))
        step_stopMin <- min(as.integer(as.character(loh_df$window_stop)))
        step_stopMax <- max(as.integer(as.character(loh_df$window_stop)))
        
        ## Calculate number of step-sized windows needed 
        intervals <- ((step_startMax + step) - step_startMin)/step
        
        step_start <- vector()
        step_stop <- vector()
        
        for (i in 1:intervals) {
            ## Obtain positions for the first step-sized window
            if (i == 1) {
                step_start[i] <- step_startMin
                step_stop[i] <- step_start[1] + step 
            }
            if (i > 1) {
                step_start[i] <- step_start[i-1] + step
                step_stop[i] <- step_stop[i-1] + step
            }
        }
        step_df <- data.frame(cbind(step_start, step_stop))
        
        ## Set the end position of the last step-sized window to be the 
        ## end of the chromosome
        step_df[nrow(step_df),2] <- as.numeric(as.character(step_stopMax))
        step_df$chromosome <- chromosome
        step_df$sample <- sample
        
        ## Get the avg loh in overlapping windows for each step
        for (w in 1:nrow(step_df)) {
            start <- as.numeric(as.character(step_df$step_start[w]))
            stop1 <- as.numeric(as.character(step_df$step_stop[w]))
            
            ## Obtain dataset with loh call positions are within the 
            ## step boundaries
            df <- loh_df[loh_df$position >= start & 
                             loh_df$position <= stop1,]
            if (nrow(df) == 0) {
                step_df$loh_avg[w] <- 0
            }
            if (nrow(df) >= 1) {
                step_df$loh_avg[w] <- mean(df$loh_diff_avg)
            }
        }
        
        ## Remove rows without loh calls made in the step-sized window
        final_step_df <- step_df[step_df$loh_avg > 0,]
        return(final_step_df)  
    })
    return(step_final)
}

#' Obtain window information
#' 
#' Calculate window positions to perform LOH calculation
#' @name lohSpec_windowPosition
#' @param out object of class dataframe with columns 'chromosome', 
#' 'position', 'n_vaf', 't_vaf', and 'sample'
#' @param step integer with the length of divisions (bp) in chromosomes
#' @param window_size integer with the size of the sliding window (bp) to be 
#' applied
#' @return list containing window start/stop positions for each chromosome 
#' from each sample to perform LOH calculations


lohSpec_windowPosition <- function(out, step, window_size)
{
    window <- lapply(out, function(x) {
        ## Calculate the number of windows necessary for each list
        min <- integer()
        max <- integer()
        window_stop_1 <- integer()
        window_num <- integer()
        min <- as.integer(min(as.numeric(as.character(x$position))))
        max <- as.integer(max(as.numeric(as.character(x$position))))
        window_stop_1 <- min+window_size
        num <- as.integer((max-min)/step)
        num <- as.vector(1:num)
        window_data_start <- vector()
        window_data_stop <- vector()
        
        ## Calculate exact window positions
        window_data <- lapply(num, function(x){
            window_data_start[x] <- as.integer(min+(step*(x-1)))
            window_data_stop[x] <- as.integer(window_stop_1+(step*(x-1)))
            window_data <- cbind(window_data_start[x], window_data_stop[x])
            return(window_data)
        })
        window_data <- plyr::ldply(window_data, data.frame)
        # Get window positions whose values are below max & set max as the 
        # final window position (end of the chromosome)
        colnames(window_data) <- c("window_start", "window_stop")
        window_final <- window_data[window_data$window_stop <= max,]
        window_final[nrow(window_final), 2] <- max
        return(window_final)
    })
    return(window)
}


