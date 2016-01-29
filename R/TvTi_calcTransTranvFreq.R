#' Calculate Transition/Transversion Frequency
#'
#' Given a data frame with columns reference, variant, sample, and trans/tranv
#' calculate the frequencies of transitions and transversion occuring.
#' @name TvTi_calcTransTranvFreq
#' @param x Object of class data frame containing columns 'reference',
#' 'variant', 'sample', 'trans_tranv'
#' @return Object of class data frame with Frequency and Proportion of
#' Transistions/Transversions appended on a sample level

TvTi_calcTransTranvFreq <- function(x)
{
    # Ensure all possible combinations of trans/tranv are represented
    trans_tranv <- c("A->C or T->G (TV)", "A->G or T->C (TI)",
                     "A->T or T->A (TV)", "G->A or C->T (TI)",
                     "G->C or C->G (TV)", "G->T or C->A (TV)")
    sample <- c('dummy_sample')
    reference <- c('A')
    variant <- c('T')
    dummy_data <- data.frame(reference, variant, sample, trans_tranv)
    x <- rbind(dummy_data, x)

    # calculate the frequency of transitions/transversions on a sample basis
    x_freq <-  table(x$trans_tranv, x$sample)

    # calculate the proportion of transitions/transversions on a sample basis
    x_prop <-  prop.table(x_freq, 2)

    # format and remove the dummy data introduced above
    x_freq <- as.data.frame(x_freq)
    x_prop <- as.data.frame(x_prop)
    x <- cbind(x_freq, x_prop$Freq)
    colnames(x) <- c('trans_tranv', 'sample', 'Freq', 'Prop')
    x <- x[which(x$sample != "dummy_sample"),]

    return(x)
}
