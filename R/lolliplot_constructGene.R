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
    
    # Construct basic gene information, if there are no domains return the gene
    gene <- data.frame(Domain=gene, pos_from=1, pos_to=length, nest=1)
    if(nrow(na.omit(domain_data)) == 0)
    {
        gene$height_min <- .1/(as.numeric(gene$nest))
        gene$height_max <- -.1/(as.numeric(gene$nest))
        gene$pos_from <- as.numeric(gene$pos_from)
        gene$pos_to <- as.numeric(gene$pos_to)
        return(gene)
    }
    
    # rename columns for domain_data and make sure description column is
    # not a factor
    colnames(domain_data) <- c("description", "start", "end")
    domain_data$description <- as.character(domain_data$description)
    
    # quality check of domain data
    if(max(domain_data$end) > length)
    {
        memo <- paste0("The end position of a domain: ",  max(domain_data$end),
                       " is exceeding the length of the protein:", length)
        warning(memo)
    } else if(min(domain_data$start) < 1) {
        memo <- paste0("The start position of a domain:",
                       min(domain_data$start),
                       "is less than the start of the protein", 1)
        warning(memo)
    }
    
    # Check that start coordinates are always less than the end coordinates
    if(any(domain_data$start >= domain_data$end))
    {
        memo <- paste0("Found a start position greater than an end position",
                       " in the protein features track. Check input to Z or",
                       "results of the biomaRt query using dataOut==TRUE.")
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
    colnames(domain_data) <- c("Domain", "pos_from", "pos_to", "nest")
    
    # combine gene and domain information
    gene <- rbind(gene, domain_data)
    
    # annotate display heights based on nesting and make sure coord are numeric
    gene$height_min <- .1/(as.numeric(gene$nest))
    gene$height_max <- -.1/(as.numeric(gene$nest))
    gene$pos_from <- as.numeric(gene$pos_from)
    gene$pos_to <- as.numeric(gene$pos_to)
    
    return(gene)
}