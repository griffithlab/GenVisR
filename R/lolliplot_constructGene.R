#' Construct gene information
#' 
#' Build gene for input into lolliplot_buildMain
#' @name lolliplot_constructGene
#' @param gene character string specifying gene name
#' @param domain_data object of class data frame specifying protien domain
#' information, obtained from lolliplot_fetchDomain, should contain columns
#' giving "description", "start", "end"
#' @param length integer specifying length of transcript in amino acids
#' @return object of class data frame giving gene and domain information

lolliplot_constructGene <- function(gene, domain_data, length)
{
    # message
    message("Constructing gene track")
    
    # rename columns for domain_data
    colnames(domain_data) <- c("description", "start", "end")
    
    # quality check of domain data
    if(max(domain_data$end) > length)
    {
        memo <- paste0("The end position of a domain:",  max(domain_data$end),
                       "is exceeding the length of the protein:", length)
        warning(memo)
    } else if(min(domain_data$start) < 1) {
        memo <- paste0("The start position of a domain:",
                       min(domain_data$start),
                       "is less than the start of the protein", 1)
        warning(memo)
    }
    
    # determine which regions are overlapping and annotate which nest domain is
    # sort on start
    domain_data$start <- as.numeric(domain_data$start)
    domain_data$end <- as.numeric(domain_data$end)
    domain_data <- domain_data[order(domain_data$start),]
    
    # annotate nests
    nest <- vector('numeric')
    end <- vector('numeric')
    for(i in 1:nrow(domain_data))
    {
        # Remove from end any values <= gene$start[i]
        idx <- domain_data$start[i] < end
        end <- end[idx]
        
        nest <- c(nest, length(end))
        end <- c(end, domain_data$end[i])
    }
    
    # add this nest information to the data frame
    domain_data$nest <- nest + 1
    
    # add in the actual transcript track to the domain information
    gene <- c(gene, 1, length, 1)
    
    # combine gene and domain information
    gene <- rbind(gene, domain_data)
    colnames(gene) <- c("Domain", "pos_from", "pos_to", "nest")
    
    # annotate display heights based on nesting and make sure coord are numeric
    gene$height_min <- .5/(as.numeric(gene$nest))
    gene$height_max <- -.5/(as.numeric(gene$nest))
    gene$pos_from <- as.numeric(gene$pos_from)
    gene$pos_to <- as.numeric(gene$pos_to)
    
    return(gene)
}