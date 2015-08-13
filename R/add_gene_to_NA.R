#' Assign NA samples gene
#' 
#' Replace NA values in a gene column to the top gene name
#' @name add_gene_to_NA
#' @param data_frame a data frame in MAF format
#' @return a data frame with NA values in a gene column coerced to the top gene
#' name

add_gene_to_NA <- function(data_frame)
{
  ##############################################################################
  # Function to replace na values in a dataframe with a gene name, (i.e.
  # when the samples with NA values for gene/trv_type are read in add a
  # gene to these NA values so that they are plotted on the x axis 
  # instead of an "NA" gene appearing)
  ##############################################################################
  
  # Get The gene with the most Mutations and add the NA samples to that gene
  # (ensures that the NA's are added in as gene with most mutations will always
  # be plotted)
  top_gene <- na.omit(rev(data_frame$gene))[1]
  data_frame$gene <- replace(data_frame$gene, is.na(data_frame$gene), top_gene)
  
  return(data_frame)
}