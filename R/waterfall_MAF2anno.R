#' Convert MAF File
#' 
#' Convert columns of a mutation annotation file "MAF" into a format
#' recognizable by internal functions
#' @name waterfall_MAF2anno
#' @param x a data frame in MAF format
#' @param label_col Character string specifying the column name of a
#' label column
#' @return a data frame coerced from MAF to TGI format

waterfall_MAF2anno <- function(x, label_col)
{
    # Check that correct column names are present and convert to internal format
    expec_col <- c('Tumor_Sample_Barcode', 'Hugo_Symbol',
                   'Variant_Classification')
    
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
    
    x <- x[,c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification',
              label_col)]
    
    if(!is.null(label_col))
    {
        colnames(x) <- c('sample', 'gene', 'trv_type', 'label')
    } else {
        colnames(x) <- c('sample', 'gene', 'trv_type')
    }
    return(x)
}