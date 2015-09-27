#' Convert MGI File
#' 
#' Convert columns of a mutation annotation file "MGI" into a format
#' recognizable by internal functions
#' @name waterfall_MGI2anno
#' @param x a data frame in MGI internal format
#' @param label_col Character string specifying the column name of a label
#' column
#' @return a data frame coerced from MGI to internal annotation format

waterfall_MGI2anno <- function(x, label_col)
{
    # Check that correct column names are present and convert to internal format
    expec_col <- c('sample', 'gene_name', 'trv_type')
    if(!is.null(label_col))
    {
        expec_col <- c(expec_col, label_col)
    }
    
    if(!all(expec_col %in% colnames(x)))
    {
        memo <- paste0("Did not detect correct column names, column names
                       should be: ", toString(expec_col))
        stop(memo)
    }
    
    x <- x[,c('sample', 'gene_name', 'trv_type', label_col)]
    if(!is.null(label_col))
    {
        colnames(x) <- c('sample', 'gene', 'trv_type', 'label')
    } else {
        colnames(x) <- c('sample', 'gene', 'trv_type')
    }
    
    return(x)
}