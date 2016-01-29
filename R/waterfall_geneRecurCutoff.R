#' Mutation Recurrence Cutoff
#'
#' Subset a MAF file keeping only samples that meet a mutation recurrence cutoff
#' @name waterfall_geneRecurCutoff
#' @param x data frame in long format with columns 'gene', 'trv_type', 'sample'
#' @param recurrence_cutoff integer specifying removal of entries not seen
#' in at least "x" percent of samples
#' @return a subset data frame
#' @importFrom plyr count

waterfall_geneRecurCutoff <- function(x, recurrence_cutoff)
{
    if(recurrence_cutoff != 0)
    {
        message("Performing recurrence cutoff...")
    }
    mutRecur <- plyr::count(unique(x), vars=c("gene"))
    mutRecur <- na.omit(mutRecur)
    mutRecur$prop <- mutRecur$freq/nlevels(x$sample)

    # If recurrence cutoff specified exceeds upper limit such that no plot
    # useful would be generated, reset recurrence cutoff
    maxRecur <- max(mutRecur$prop)
    if(maxRecur < recurrence_cutoff)
    {
        memo <- paste0("The recurrence cutoff specified exceeds the recurrence",
                       " seen in the data, resetting this value to equal max ",
                       "recurrence:", maxRecur)
        warning(memo)
        recurrence_cutoff <- maxRecur
    }

    gene_above_recur <- mutRecur[mutRecur$prop >= recurrence_cutoff,]$gene
    # add NA to the end of 'gene_above_recurrence' vector, allowing for all
    # samples having NA as a gene name to be retained in the subset below
    gene_above_recur <- c(as.character(gene_above_recur), NA)

    # subset the original data frame based on the following: keep gene if it is
    # in the gene vector in "mutation_recurrence_subset"
    x <- x[(x$gene %in% gene_above_recur), ]

    return(x)
}
