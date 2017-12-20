################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#' Class LOH
#' 
#' An S4 class for the lohSpec plot object
#' @name lohSpec-class
#' @rdname lohSpec-class
#' @slot lohFreq_plot gtable object for the lohFreq plot
#' @slot lohSpec_plot gtable object for the lohSpec plot
#' @slot lohData data.table object soring loh data with column names: sample, 
#' chromosome, position, t_vaf, n_vaf. 
#' @exportClass lohSpec
#' @importFrom data.table data.table
#' @importFrom gtable gtable
methods::setOldClass("gtable")
setClass(
    Class="lohSpec",
    representation=representation(lohFreq_plot="gtable",
                                  lohSpec_plot="gtable",
                                  lohData="data.table"),
    validity = function(object) {
        
    }
)

#' Constructor for the lohSpec class
#' 
#' @name lohSpec
#' @rdname lohSpec-class
#' @param input Object of class VarScan.
#' @param Character vector specifying the chromosomes of interest. If NULL,
#' will use autosomes for human (chr1-22). 
#' @param samples Character vector specifying samples to plot. If not NULL
#' all samples in "input" not specified with this parameter are removed.
#' @param boundaries Object of class data frame with rows representing chromosome
#' boundaries for a genome assembly. The data frame must contain columns with
#' the following names "chromosome", "start", "end". If let null, will determine
#' chr boundaries using preloaded/specified genome.
#' @param genome Character string specifying a valid UCSC genome (see details).
#' @param gender Character vector of length equal to the number of samples,
#' consisting of elements from the set {"M", "F"}. Used to suppress the plotting
#' of allosomes where appropriate.
#' @param step Integer value specifying the step size (i.e. the number of base
#' pairs to move the window). required when method is set to slide
#' (see details).
#' @param window_size Integer value specifying the size of the window in base
#' pairs in which to calculate the mean Loss of Heterozygosity (see details).
#' @param normal Numeric value within the range 0-1 specifying the expected
#' normal variant allele frequency to be used in Loss of Heterozygosity 
#' calculations. defaults to .50\%


lohSpec <- function(input, chr=NULL, samples=NULL, y=NULL, genome='hg19',
                    gender=NULL, step=1000000, window_size=2500000, 
                    normal=.50, gradient_midpoint=.2, gradient_low="#ffffff",
                    gradient_mid="#b2b2ff", gradient_high="#000000",
                    theme_layer=NULL,  verbose){
    ## Calculate all data for plots
    loh_data <- lohData(input, chr=chr, samples=samples, y=y, genome=genome,
                    step=step, window_size=window_size,
                    normal=normal, verbose)
    
    ## Use the lohData to generate lohSpec plots
    lohSpec_plot <- lohSpec_buildMainPlot(object=loh_data, plotLayer=NULL)
    
    
}
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Private Classes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Private Class lohData
#' 
#' An S4 class for the Data of the loh plot object
#' @name lohData-class
#' @name lohData-class
setClass("lohData",
         representation=representation(primaryData="data.table",
                                       windowData="data.table",
                                       windowCalcData="data.table",
                                       chrData="data.table"),
         validity = function(object){
             
         }
         )

#' Constructor for the lohData class.
#' 
#' @name lohData
#' @rdname lohData-class
#' @param object Object of class VarScan 
lohData <- function(object, chr, samples, y, genome, step, window_size, 
                    normal, verbose) {
    ## Get the primary loh data
    object <- VarScanFormat(path = "~/Google Drive/HCC1395.varscan.tsv")
    primaryData <- getLohData(object=object, chr=chr, verbose = verbose)
    
    ## Quality check on the primary data
    ## To-Do: Put the quality check in the validity check
    primaryData <- lohSpec_qual(object=primaryData)
    
    library(BiocInstaller)
    biocLite("BSgenome.Hsapiens.UCSC.hg19")
    BSgenome <- getBSgenome(genome = "BSgenome.Hsapiens.UCSC.hg19")
    
    ## Get the chromosome data
    if(is.null(y)) {
        y <- data.table()
    }
    ## Check that y is a data.table
    if (!is.data.table(y)) {
        message("y is not a data.table, attempting to coerce")
        y <- data.table(y)
    }
    chrData <- getChrBoundaries(object=y, genome=genome)
   
    ## Produce data.table with window position data
    window_data <- getLohSlidingWindow(object = primaryData, step = step, 
                                       window_size = window_size)
    
    ## Perform loh calculations on each chromosome and sample within each window
    loh_abs_diff <- getLohCalculation(object=primaryData, 
                                      window_data=window_data, normal=normal)
    
    ## Calculate avg loh for overlapping regions
    loh_abs_diff_overlap <- rbindlist(getLohStepCalculation(object=loh_abs_diff, step=step))
   
    ## Initialize the object
    new("lohData", primaryData=primaryData, windowData=rbindlist(window_data), 
        windowCalcData=loh_abs_diff_overlap, chrData=chrData)
}


################################################################################
###################### Accessor function definitions ###########################

#########################################################
##### Function to perform quality check on loh data #####
#' @rdname lohSpec_qual-methods
#' @param object of class lohData
#' @aliases lohSpec_qual
#' @return Data.table with quality check 
setMethod(f="lohSpec_qual",
          signature="data.table", 
          definition=function(object){
              primaryData <- object
              ## Check that values supplied in vaf columns are in the expected range
              if (any(primaryData$tumor_var_freq>1 | primaryData$normal_var_freq>1)) {
                  stop("Detected values in either the normal or tumor variant ",
                       "allele fraction columns above 1. Values supplied should ",
                       "be a proportion between 0-1!")
              }
              if (any(primaryData$normal_var_freq<0.4 | primaryData$normal_var_freq>0.6)){
                  message("Detected values with a variant allele fraction either ",
                          "above .6 or below .4 in the normal. Please ensure ",
                          "variants supplied are heterozygous in the normal!")
                  message("Removing coordinates with normal VAF > 0.6 or < 0.4")
                  primaryData <- primaryData[normal_var_freq<=0.6 & 
                                                 normal_var_freq>=0.4]
              }
              
              ## Check the chromosome column - see if it has "chr" as the prefix
              if(!all(grepl("^chr", primaryData$chrom))) {
                  memo <- paste0("Did not detect the prefix chr in the chromosome column",
                                 " of x... adding prefix")
                  message(memo)
                  primaryData$chrom <- paste0("chr", primaryData$chrom)
              } else if (all(grepl("^chr", primaryData$chrom))) {
                  message(paste0("detected chr in the chromosome column of x...",
                                 "proceeding"))
              } else {
                  stop("Detected unknown or mixed prefixes in the chromosome",
                       " column of x... should either be chr or none i.e. ",
                       "chr1 or 1")
              }
              
              ## Change column names
              colnames(primaryData) <- c("chromosome", "position", "t_vaf", 
                                         "n_vaf", "sample")
              return(primaryData)
          })

##### FIX THIS TO USE DATA.TABLES #####
#####################################################
##### Function to get the chromosome boundaries #####
#' @rdname getChrBoundaries-methods
#' @param object of class lohData 
#' @param genome character specifying which genome to use
#' @return Data.table with chr and start/stop positions
#' @aliases getChrBoundaries
setMethod(f="getChrBoundaries", 
          signature="data.table",
          definition=function(object, genome){
              ## Preloaded genome options
              preloaded <- c('hg38', 'hg19', 'mm10', 'mm9', 'rn5')
              if (nrow(object) == 0 ) {
                  ## Check that genome specified is not the ensembl name
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
                      message("attempting to query UCSC sql database for chromosome
                          positions")
                      cyto_data <- suppressWarnings(multi_cytobandRet(genome))
                      chr_pos <- multi_chrBound(cyto_data)
                  }
                  if (any(genome == preloaded)) {
                      message("genome specified is preloaded, retrieving data...")
                      chr_pos <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == genome,]
                      chr_pos <- multi_chrBound(chr_pos)  
                  }
              } 
              if (nrow(object) != 0){
                  if(!all(c('chromosome', 'start', 'end') %in% colnames(y)))
                  {
                      memo <- paste0("Did not detect correct column names in y, missing",
                                     "one of \"chromosome\", \"start\", \"end\"")
                      stop(memo)
                  }
                  # Ensure that columns in data frame are of proper type
                  object$chromosome <- as.character(object$chromosome)
                  object$start <- as.integer(as.character(object$start))
                  object$end <- as.integer(as.character(object$end))
                  message("detected input to y, using supplied positions for chromosome
                          boundaries")
                  chr_pos <- object
              }
              
              # Quality check for dummy data
              if(nrow(chr_pos) < 1)
              {
                  memo <- paste0("Could not retrieve chromosome boundaries from",
                                 " UCSC, please specify this information via ",
                                 "the y paramter")
                  stop(memo)
              }
              return(data.table(chr_pos))
          })

##########################################################################
##### Function to generate window position data for loh calculations #####
#' @rdname getLohSlidingWindow-methods
#' @param object of class lohData 
#' @param step integer specifying the step size between the start position of
#' each window
#' @param window_size integer specifying the window size for loh calcuations
#' @return Data.table with window start/stop positions
#' @aliases getLohSlidingWindow
setMethod(f="getLohSlidingWindow",
          signature="data.table",
          definition=function(object, step, window_size, ...){
              object <- primaryData
              ## Obtain lists for each sample and chromosome
              out <- split(object, list(as.character(object$chromosome),
                                          as.character(object$sample)))
              
              ## Obtain the window position values
              window <- lapply(out, function(x, step, window_size) {
                  ## Get the min and max position on the chromosome
                  min <- integer()
                  max <- integer()
                  window_stop_1 <- integer()
                  window_num <- integer()
                  min <- as.integer(min(as.numeric(as.character(x$position))))
                  max <- as.integer(max(as.numeric(as.character(x$position))))
                  ## Get the end of the first window position
                  window_stop_1 <- min+window_size
                  ## Calculate the number of windows necessary
                  num <- as.integer((max-min)/step)
                  num <- as.vector(1:num)
                  window_data_start <- vector()
                  window_data_stop <- vector()
                  
                  ## Calculate exact window positions
                  window_data <- lapply(num, function(x){
                      window_data_start[x] <- as.integer(min+(step*(x-1)))
                      window_data_stop[x] <- as.integer(window_stop_1+(step*(x-1)))
                      window_data <- data.table(cbind(window_data_start[x], window_data_stop[x]))
                      return(window_data)
                  })
                  window_data <- rbindlist(window_data)
                  # Get window positions whose values are below max & set max as the 
                  # final window position (end of the chromosome)
                  colnames(window_data) <- c("window_start", "window_stop")
                  window_final <- window_data[window_data$window_stop <= max,]
                  window_final[nrow(window_final), 2] <- max
                  ## Put in the chromosome 
                  window_final$chromosome <- as.character(x$chromosome[1])
                  return(window_final)
              }, 
              step = step, window_size = window_size)
              
              return(window)
          })

###############################################################
##### Function to perform loh calcluations in each window #####
#' @rdname getLohCalculation-methods
#' @param object of class lohData 
#' @param window_data of class data.table 
#' @param normal integer specifying normal vaf
#' @aliases getLohCalculation
setMethod(f="getLohCalculation", 
          signature="data.table",
          definition=function(object, window_data, normal, ...) {
              object <- split(object, list(as.character(object$chromosome),
                                           as.character(object$sample)))
              window_data <- window_data
              ## Separate out sample and window data by chromosome name
              df <- lapply(object, function(sample_data, window, 
                                                      normal) {
                  chromosome <- as.character(sample_data[1,1])
                  sample <- as.character(sample_data[1,5])
                  chromosome.sample <- paste(chromosome, sample, sep = ".")
                  window <- window_data[[grep(chromosome.sample, names(window_data))]]
                  ## For each window position, get the vaf data that falls 
                  ## within that window
                  dataset <- rbindlist(apply(window, 1, function(x, 
                                                                 sample_data, normal){
                      if (x[3] != as.character(sample_data[1,1])) {
                          stop("Chromosomes in window and sample vaf data do not match")
                      }
                      w_start <- as.numeric(as.character(x[1]))
                      w_stop <- as.numeric(as.character(x[2]))
                      ## Filter out vaf data outside the window
                      filtered_data <- sample_data[position >= w_start &
                                                       position <= w_stop]
                      
                      ## Peroform loh calclulation to obtain avg loh in the 
                      ## window's frame
                      loh_calc_avg <- mean(abs(as.numeric(as.character(
                          filtered_data$t_vaf)) - normal))
                      if (is.na(loh_calc_avg)) {
                          loh_calc_avg <- NA
                          w_start <- NA
                          w_stop <- NA
                      }
                      filtered_data$loh_diff_avg <- loh_calc_avg
                      filtered_data$window_start <- w_start
                      filtered_data$window_stop <- w_stop
                      return(filtered_data)
                  }, 
                  sample_data=sample_data, normal=normal))
                  dataset <- na.omit(dataset, cols = c("loh_diff_avg", 
                                                       "window_start", 
                                                       "window_stop"))
                  return(dataset)
              }, window=window_data, normal=normal)
              return(df)
          })

#######################################################################
##### Function to perform loh calcluations in overlapping windows #####
#' @rdname getLohStepCalculation-methods
#' @param object of class lohData
#' @param step integer 
#' @aliases getLohStepCalculation
setMethod(f = "getLohStepCalculation",
          signature="list",
          definition=function(object, step, ...) {
              object <- loh_abs_diff
              step_loh_calc <- lapply(object, function(x, step) {
                  ## Get the sample and chromosome information
                  sample <- unique(x$sample)
                  chromosome <- unique(x$chromosome)
                  
                  ## Obtain boundaries for each step-sized window
                  start <- unique(x$window_start)
                  stop <- c(start[-1], max(x$window_stop))
                  step_boundaries <- data.table(chromosome=chromosome, start=start, stop=stop)
                  step_boundaries$sample <- sample
                  
                  ## Get the average loh within each step-sized window
                  loh_df <- x
                  loh_step_avg <- apply(step_boundaries, 1, function(x, loh_df_data) {
                      start <- as.numeric(as.character(x[2]))
                      stop <- as.numeric(as.character(x[3]))
                      step_df <- loh_df_data[position >= start & 
                                      position < stop]
                      if (nrow(step_df) == 0) {
                          loh_step_avg <- 0
                          }
                      if (nrow(step_df) > 0) {
                          loh_step_avg <- mean(step_df$loh_diff_avg)
                          }
                      return(loh_step_avg)
                  }, loh_df_data=loh_df)
                  step_boundaries$loh_step_avg <- loh_step_avg
                  return(step_boundaries)
              }, step=step)
              return(step_loh_calc)
          })

#######################################################################
##### Function to perform loh calcluations in overlapping windows #####
#' @rdname lohSpec_buildMainPlot-methods
#' @param object of class lohData
#' @param step integer 
#' @aliases lohSpec_buildMainPlot
setMethod(f = "lohSpec_buildMainPlot",
          signature="lohData",
          definition=function(object, ...) {
              x <- object@windowCalcData
              x <- x[loh_step_avg > 0]
              
              ## Set the order of the chromosomes
              chr <- gtools::mixedsort((unique(x$chromosome)))
              sample <- gtools::mixedsort((unique(x$sample)))
              x$chromosome_f <- factor(x$chromosome, levels=chr)
              x$sample <- factor(x$sample, levels=sample, labels=sample)
              
              dummyData <- object@chrData
              # define dummy data which will be chromosome boundaries, these are plotted
              # but are transparent and will not appear in the plot
              dummy_data <- geom_rect(data=dummyData, aes_string(xmin='start', xmax='end',
                                                                 ymin=-1, ymax=1),alpha=0)
              
              # Define the main plot
              data <- geom_rect(data=x, aes_string(xmin='start',
                                                   xmax='stop',
                                                   ymin=-1,
                                                   ymax=1, fill='loh_step_avg'))
              
              # Define additional plot parameters
              facet <- facet_grid(sample ~ chromosome_f, scales="free", space="free")
              
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
              
              LOHgradient <- scale_fill_gradient2(midpoint = gradient_midpoint,
                                                  guide="colourbar",
                                                  high=gradient_high,
                                                  mid=gradient_mid,
                                                  low=gradient_low,
                                                  space='Lab')
              
              # Build the plot
              tmp <- data.frame(x=0, y=0)
              p1 <- ggplot(data=tmp, aes(y=0)) + dummy_data + data + facet + x_scale + y_scale + 
                  lab_x + lab_y + BWscheme + LOHgradient + plotTheme + plotLayer
              print(p1)
              return(p1)
          })

