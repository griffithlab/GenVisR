#' Convert Custom File
#' 
#' Convert columns of a Custom annotation file into a format
#' recognizable by internal functions
#' @name waterfall_Custom2anno
#' @param x a data frame with columns having values for sample, gene, mutation
#' type
#' @param label_col Character string specifying the column name of a
#' label column (optional)
#' @return a data frame coerced from custom to annotation format

waterfall_Custom2anno <- function(x, label_col)
{
    # message statement
    memo <- paste0("Detected \"Custom\" file_type flag, ",
                   "looking for correct column names...")
    message(memo)
    
    # define expected columns
    expec_col <- c("sample", "gene", "variant_class")
    if(!is.null(label_col))
    {
        expec_col <- c(expec_col, label_col)
    }
    
    # check expected columns are present
    if(!all(expec_col %in% colnames(x)))
    {
        stop("Did not detect correct column names, check file_type flag?")
    }
    
    x <- x[,c('sample', 'gene', 'variant_class', label_col)]
    
    if(!is.null(label_col))
    {
        colnames(x) <- c('sample', 'gene', 'trv_type', 'label')
    } else {
        colnames(x) <- c('sample', 'gene', 'trv_type')
    }
    
    # if no silent mutations are present warn the user
    if(all(!toupper(x$trv_type) %in% toupper("silent")))
    {
        warning("Did not detect silent mutations in input, is this expected?")
    }
    return(x)
}