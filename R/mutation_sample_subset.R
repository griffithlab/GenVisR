#' Mutation Sample Cutoff
#' 
#' Subset a MAF file keeping only entries in a selection of genes
#' @name mutation_sample_subset
#' @param x a data frame in long format with columns 'gene', 'trv_type'
#' @param genes character vector listing genes to plot
#' @return a subset data frame

mutation_sample_subset <- function(x, genes)
{
  # Perform quality checks
  if(typeof(genes) != 'character' & class(genes) != 'character')
  {
    warning("argument supplied to main.genes is not a character vector, attempting to coerce")
    genes <- as.character(genes)
  }
  if(!all(toupper(genes) %in% toupper(x$gene)))
  {
    warning("genes supplied in main.genes contains an element not found in x")
  }
  
  genes <- c(genes, NA)
  x <- x[(toupper(x$gene) %in% toupper(genes)), ]
  
  return(x)
}