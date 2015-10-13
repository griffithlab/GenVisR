#' Check input to lohView
#'
#' Perform data quality checks on input supplied to lohView
#' @name lohView_qual
#' @param x object of class data frame with columns 'chromosome', 'position',
#' 'n_vaf', 't_vaf', 'sample'
#' @return list of inputs passing basic quality controls

lohView_qual <- function(x)
{
    # check input data to x
    if(!is.data.frame(x))
    {
        stop("Did not detect a data frame for input to x")
    }

    # check that correct columns are supplied in x
    x_col <- c('chromosome', 'position', 'n_vaf', 't_vaf', 'sample')
    if(!all(x_col %in% colnames(x)))
    {
        stop('Did not detect required column names in x, required columns are: '
        , paste0(x_col, sep="\t"))
    }

    return(list(x))
}
