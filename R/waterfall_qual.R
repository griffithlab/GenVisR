#' Check input to mutSpec
#'
#' Perform a data quality check on input to mutSpec
#' @name waterfall_qual
#' @param x a data frame in annotation format
#' @param y a data frame containing clinical data or a null object
#' @param z a data frame containing mutation burden information or a null object
#' @param file_type Character string specifying the input format to expect in x
#' @param label_col Character string specifying the column name of a label
#' @param variant_class_order Character vector specifying the hierarchical order
#' of mutation types to plot
#' @return a list of data frames passing quality checks

waterfall_qual <- function(x, y, z, file_type, label_col, variant_class_order)
{
    # print message statement
    message("Checking if input is properly formatted...")

    # Check input data to x
    if(!is.data.frame(x))
    {
        stop("Did not detect a data frame for input to x")
    }

    # Convert file type to internal format
    if(toupper(file_type) == toupper("MAF"))
    {
        x <- waterfall_MAF2anno(x, label_col, variant_class_order)
    } else if(toupper(file_type) == toupper("MGI")) {
        x <- waterfall_MGI2anno(x, label_col, variant_class_order)
    } else if(toupper(file_type) == toupper("Custom")) {
        x <- waterfall_Custom2anno(x, label_col, variant_class_order)
    } else if(toupper(file_type) == toupper("VEP")) {
        x <- waterfall_VEP2anno(x, label_col, variant_class_order)
    } else {
        stop("Unrecognized file_type: ", file_type)
    }
    
    # drop unused levels in x
    x$sample <- as.factor(x$sample)
    x$gene <- as.factor(x$gene)
    x$trv_type <- as.factor(x$trv_type)
    x <- droplevels(x)

    # Check input data to clinDat
    if(!is.null(y))
    {
        if(!is.data.frame(y))
        {
            stop("Did not detect a data frame for input to clinDat")
        }
        y <- droplevels(y)

        if(!all(c('sample', 'variable', 'value') %in% colnames(y)))
        {
            stop("Did not detect correct sample names in clinDat")
        }
        
        # make sure clinical data columns are of expected class
        y$sample <- as.factor(y$sample)
        y$variable <- as.factor(y$variable)
        y$value <- as.character(y$value)
    }

    # check input data to mutBurden
    if(!is.null(z))
    {
        if(!is.data.frame(z))
        {
            stop("Did not detect a data frame for input to mutBurden")
        }
        z <- droplevels(z)

        if(!all(c('sample', 'mut_burden') %in% colnames(z)))
        {
            stop("Did not detect correct sample names in mutBurden")
        }
        
        # Make sure mutation burden columns are of proper class
        z$sample <- as.factor(z$sample)
        z$mut_burden <- as.numeric(as.character(z$mut_burden))
    }
    
    return(list(x, y, z))
}
