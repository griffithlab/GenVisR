#' Convert Custom File
#'
#' Convert columns of a Custom annotation file into a format
#' recognizable by internal functions
#' @name waterfall_Custom2anno
#' @param x Object of class data frame with rows representing mutations and
#' containing columns "key, "sample", "gene", "variant_class".
#' @param label_col Character string specifying the column name of a
#' label column (optional)
#' @return a data frame coerced from custom to annotation format

waterfall_Custom2anno <- function(x, label_col, variant_class_order)
{
    # message statement
    memo <- paste0("Detected \"Custom\" fileType flag, ",
                   "looking for correct column names...")
    message(memo)

    # define expected columns
    expec_col <- c("key", "sample", "gene", "variant_class")
    if(!is.null(label_col))
    {
        expec_col <- c(expec_col, label_col)
    }

    # check expected columns are present
    if(!all(expec_col %in% colnames(x)))
    {
        memo <- paste0("Did not detect correct column names, column names
                       should be: ", toString(expec_col))
        stop(memo)
    }

    x <- x[,c("key", 'sample', 'gene', 'variant_class', label_col)]
    if(!is.null(label_col))
    {
        colnames(x) <- c('sample', 'gene', 'trv_type', 'label')
    } else {
        colnames(x) <- c('sample', 'gene', 'trv_type')
    }
    
    if(is.null(variant_class_order)){
        memo <- paste0("Detected NULL in variant_class_order, ",
                       "this parameter is required if fileType is set ",
                       "to \"Custom\"")
        stop(memo)
    }
    
    # Check that elements in trv_type are in the mutation order
    if(any(!x$trv_type %in% mutation_order))
    {
        memo <- paste0("Detected an invalid mutation type, valid values for ",
                       "Custom are: ", toString(mutation_order))
        stop(memo)
    }
    
    # remove duplicate keys keeping entries based on a trv_type hiearchy
    x$trv_type <- factor(x$trv_type, levels=variant_class_order)
    x <- x[order(x$sample, x$key, x$gene),]
    x <- x[!duplicated(x[, c("sample", "key", "gene")]),]
    
    return(x)
}
