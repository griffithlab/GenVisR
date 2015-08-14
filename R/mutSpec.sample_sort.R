#' sort samples in a MAF file
#' 
#' perform a hiearchial sort on samples based on the presence of mutations in
#' an ordered list of genes
#' @name mutSpec.sample_sort
#' @param x a data frame in long format with column names "sample", 
#' "gene", "trv_type"
#' @return a vector of samples in a sorted order

mutSpec.sample_sort <- function(x)
{
    # recast the data going from long format to wide format, values in this data are counts of a mutation call
    wide_data <- reshape2::dcast(x, sample ~ gene, fun.aggregate = length,
                                 value.var="trv_type")
    
    # apply a boolean function to convert the data frame values to 1's and 0's
    wide_data[,-1] <- t(apply(wide_data[,-1], 1, mutSpec.convert_to_boolean))
    wide_boolean <- wide_data
    
    # move the first column (sample names) to row names
    wide_boolean_format <- data.frame(wide_boolean[,-1],
                                      row.names=wide_boolean[,1])
    
    # reverse the columns so that genes with highest mutation's are listed first
    # (assuming gene_sort has been run on the data frame)
    test <- wide_boolean_format[,rev(1:ncol(wide_boolean_format))]
    
    # if there are any NA values present in a sample at the gene put that
    # NA gene last so samples with NA are plotted last
    
    if(any(colnames(test) == 'NA.'))
    {
        # Find which column has the NA header
        NA_index <- which(colnames(test) == 'NA.')
        
        # Append NA column to end of data frame
        NA_gene <- test[,NA_index]
        test <- test[,-NA_index]
        test <- cbind(test, NA_gene)
    }
    
    # hiearchial sort on all column's (i.e. genes) such that samples are
    # rearranged if there is a mutation in that gene
    
    sample_order <- rownames(test[do.call(order, as.list(-test)),])
    
    return(sample_order)
}