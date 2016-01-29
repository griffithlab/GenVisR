#' Annotate Transitions and Transversions
#'
#' Given a data frame with columns reference and variant annotate the base
#' change occurring
#' @name TvTi_annoTransTranv
#' @param x Object of class data frame containing columns 'reference', 'variant'
#' @return Object of class data frame with transition/transversion annotations
#' appended

TvTi_annoTransTranv <- function(x)
{
    # add an extra column with the reference to variant
    x$base_change <- paste0(toupper(x$reference), "2", toupper(x$variant))

    # annotate the grouping of the base change
    x$trans_tranv <- switch(x$base_change, A2C="A->C or T->G (TV)",
                            T2G="A->C or T->G (TV)", A2G="A->G or T->C (TI)",
                            T2C="A->G or T->C (TI)", A2T="A->T or T->A (TV)",
                            T2A="A->T or T->A (TV)", G2A="G->A or C->T (TI)",
                            C2T="G->A or C->T (TI)", G2C="G->C or C->G (TV)",
                            C2G="G->C or C->G (TV)", G2T="G->T or C->A (TV)",
                            C2A="G->T or C->A (TV)")

    # remove the temp base change column
    x$base_change <- NULL
    return(x)
}
