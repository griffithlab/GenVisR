#' Construct gene information
#' 
#' Build gene for input into build_lolli
#' @name lolliplot_constructGene
#' @param gene character string specifying gene name
#' @param domain_data object of class data frame specifying protien domain
#' information, obtained from lolliplot_fetchDomain, should contain columns
#' giving "description", "start", "end"
#' @param length integer specifying length of transcript in amino acids
#' @return object of class data frame giving gene and domain information

lolliplot_constructGene <- function(gene, domain_data, length)
{
    # rename columns for domain_data
    colnames(domain_data) <- c("description", "start", "end")
    
    # construct gene row for data frame
    gene <- c(gene, 1, length)
    
    # combine gene and domain information
    gene <- rbind(gene, domain_data)
    
    # add in heights for gene display and convert positions to numeric class
    gene$start <- as.numeric(gene$start)
    gene$end <- as.numeric(gene$end)
    
    # change order of positions so plot will plot largest first and determine
    # which domains are nested within other domains
    gene <- gene[order(gene$start),]
    nest <- vector('numeric')
    end <- vector('numeric')
    for(i in 1:nrow(gene))
    {
        # Remove from end any values <= gene[i]$start
        idx <- gene$start[i] < end
        end <- end[idx]
        
        nest <- c(nest, length(end))
        end <- c(end, gene$end[i])
    }
    
    gene$nest <- nest + 1
    
    gene$height_min <- .5/(gene$nest)
    gene$height_max <- -.5/(gene$nest)
    # relabel column names
    colnames(gene) <- c("Domain", "pos_from", "pos_to", "nest", "height_min",
                        "height_max")
    
    return(gene)
}