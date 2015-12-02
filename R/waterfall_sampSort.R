#' sort samples in a MAF file
#'
#' perform a hiearchial sort on samples based on the presence of mutations in
#' an ordered list of genes
#' @name waterfall_sampSort
#' @param x a data frame in long format with column names "sample",
#' "gene", "trv_type"
#' @return a vector of samples in a sorted order
#' @importFrom reshape2 dcast

waterfall_sampSort <- function(x)
{
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
    # NA gene last so samples with NA are plotted last
    if(any(grepl("^NA$", colnames(wide_boolean))))
    {
        # Find which column has the NA header
        NA_index <- which(grepl("^NA$", colnames(wide_boolean)))

        # Append NA column to end of data frame
        NA_gene <- wide_boolean[,NA_index]
        wide_boolean <- wide_boolean[,-NA_index]
        wide_boolean <- cbind(wide_boolean, NA_gene)
    }

    # hiearchial sort on all column's (i.e. genes) such that samples are
    # rearranged if there is a mutation in that gene
    sample_order <- wide_boolean[do.call(order, as.list(-wide_boolean[2:ncol(wide_boolean)])),]$sample

    # Put those samples not in sample order in from the original levels of the
    # data
    not_in <- as.character(levels(sample_order)[which(!levels(sample_order) %in% sample_order)])
    sample_order <- c(as.character(sample_order), not_in)

    return(sample_order)
}
