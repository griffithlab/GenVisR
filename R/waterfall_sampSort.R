#' sort samples in an internal waterfall file.
#'
#' perform a hiearchial sort on samples based on the presence of mutations in
#' an ordered list of genes if a sample order is unspecified.
#' @name waterfall_sampSort
#' @param x a data frame in long format with column names "sample",
#' "gene", "trv_type"
#' @param sampOrder Character vector specifying the order of samples to plot.
#' @return a vector of samples in a sorted order
#' @importFrom reshape2 dcast

waterfall_sampSort <- function(x, sampOrder=NULL)
{
    # if a sample order is already defined plot that instead
    if(!is.null(sampOrder))
    {
        sampOrder <- as.character(unique(sampOrder))
        # determine if there are any new samples and give warning if true
        new_samples <- x$sample[!x$sample %in% sampOrder]
        if(length(new_samples) != 0)
        {
            memo <- paste0("The following samples were not detected in the ",
                           "original data: ", toString(new_samples), " adding",
                           " these to the plot!")
            warning(memo)
        }
        
        # return what was given originally
        return(sampOrder)
    }
    
    # recast the data going from long format to wide format, values in this data
    # are counts of a mutation call
    wide_data <- reshape2::dcast(x, sample ~ gene, fun.aggregate = length,
                                 value.var="trv_type")

    # apply a boolean function to convert the data frame values to 1's and 0's
    values <- wide_data[,-1, drop=FALSE]
    sample <- wide_data[,1]
    values <- data.frame(apply(values, 2,
                               function(x) as.numeric(as.logical(x))))
    wide_boolean <- cbind(sample, values)

    # reverse the columns so that genes with highest mutation's are listed first
    # (assumes gene_sort has been run on the data frame)
    wide_boolean <- wide_boolean[,c(1, rev(2:ncol(wide_boolean)))]

    # if there are any NA values present in a sample at the gene put that
    # remove that sample and save (at this stage it should be samples)
    if(any(grepl("^NA.$", colnames(wide_boolean))))
    {
        # Find which column has the NA header
        NA_index <- which(grepl("^NA.$", colnames(wide_boolean)))

        # Append NA column to end of data frame
        NA_gene <- wide_boolean[,NA_index]
        wide_boolean <- wide_boolean[,-NA_index]
        
        # Save copy and remove samples with no mutations,
        # these will be added to the end
        samp_no_mut <- wide_boolean[rowSums(wide_boolean[2:ncol(wide_boolean)]) == 0,]$sample
        samp_no_mut <- as.character(samp_no_mut)
        wide_boolean <- wide_boolean[!wide_boolean$sample %in% samp_no_mut,]
    } else {
        samp_no_mut <- NULL
    }

    # hiearchial sort on all column's (i.e. genes) such that samples are
    # rearranged if there is a mutation in that gene
    sample_order <- wide_boolean[do.call(order, as.list(-wide_boolean[2:ncol(wide_boolean)])),]$sample

    # Put those samples not in sample order in from the original levels of the
    # data (these are samples with no mutations)
    not_in <- as.character(levels(sample_order)[which(!levels(sample_order) %in% sample_order)])
    not_in <- not_in[!not_in %in% samp_no_mut]
    sample_order <- c(as.character(sample_order), as.character(not_in))
    
    # Put those samples with no mutations back in
    if(!is.null(samp_no_mut))
    {
        sample_order <- c(sample_order, samp_no_mut)
    }

    return(sample_order)
}
