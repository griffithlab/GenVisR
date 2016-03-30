#' Convert VEP File
#'
#' Convert columns of a variant effect predictor file "VEP" into a format
#' recognizable by internal functions
#' @name waterfall_VEP2anno
#' @param x a data frame in VEP format
#' @param label_col Character string specifying the column name of a
#' label column
#' @return a data frame coerced from VEP to intermal expected format

waterfall_VEP2anno <- function(x, label_col)
{
    # Split fields in the "Extra" column of a VEP file into actual columns
    
    
    # Check that correct column names are present and convert to internal format
    expec_col <- c('IND', 'Gene', 'Consequence')
    
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