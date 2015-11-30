#' XY Numeric Conversion
#' 
#' Convert X and Y chromosome names to numeric 23 and 24
#' @name xy_conv
#' @param x object of class data frame with columns 'chromosome', 'position', 
#' 'n_vaf', 't_vaf', 'sample'
#' @examples
#' data <- xy_conv("HCC1_LOH.tsv")
#' @return Dataframe
#' @export

xy_conv <- function(x) {
    LOH_data <- read.delim(x, header=TRUE)
    colnames(LOH_data) <- c("chromosome", "position", "n_vaf", "t_vaf", 
                            "sample")
    LOH_data$chrom <- as.character(LOH_data$chrom)
    LOH_data$chrom[LOH_data$chrom=="X"] <- 23
    LOH_data$chrom[LOH_data$chrom=="Y"] <- 24
    LOH_data$chrom <- as.numeric(as.character(LOH_data$chrom))
    return(LOH_data)
}
