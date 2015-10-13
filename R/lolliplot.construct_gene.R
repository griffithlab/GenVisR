#' Construct gene information
#'
#' Build gene for input into build_lolli
#' @name lolliplot.construct_gene
#' @param gene character string specifying gene name
#' @param domain_data object of class data frame specifying protien domain
#' information, obtained from fetchDomain
#' @param length integer specifying length of transcript in amino acids
#' @return object of class data frame giving gene and domain information

lolliplot.construct_gene <- function(gene, domain_data, length)
{
    # construct gene row for data frame
    gene <- c(gene, 1, length)

    # combine gene and domain information
    gene <- rbind(gene, domain_data)

    # add in heights for gene display and convert positions to numeric class
    gene$height_min <- -.5
    gene$height_max <- .5
    gene$pos_from <- as.numeric(gene$pos_from)
    gene$pos_to <- as.numeric(gene$pos_to)

    # relabel column names
    colnames(gene) <- c("Domain", "pos_from", "pos_to", "height_min",
    "height_max")

    return(gene)
}
