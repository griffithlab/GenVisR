#' Construct Lolliplot from data
#' 
#' Construct a Lolliplot from object of class data frame giving observed mutations and an ensembl transcript id
#' @name lolliplot
#' @param data object of class data frame containing columns transcript_name, gene, and amino_acid_change and rows denoting mutations
#' @param cosmic boolean value specifying if cosmic mutations should be retrieved and plotted
#' @param fill_value character string giving the name of the column to shade variants on "required"
#' @param label_column character string specifying column containing text information to be plotted, defaults to NULL
#' @param plot_text_angle integer specifying angle of text to be plotted if label_column is specified
#' @param plot_text_size integer specifying the size of the text plotted if label_column is specified
#' @param point_size integer specifying the size of points plotted
#' @param gene_colour color specifying the background fill of the plotted gene
#' @param spread_degree integer specifying the amount of distance to dodge points on the x axis
#' @return object of class ggplot2
#' @export

lolliplot <- function(data, cosmic=TRUE, fill_value='trv_type', label_column=NULL, plot_text_angle=45, plot_text_size=5, point_size=1, gene_colour='#999999', spread_degree=15)
{
  # extract transcript id
  transcriptID <- as.character(data$transcript_name[1])
  
  # extract HUGO gene name
  gene <- as.character(data$gene[1])
  
  # obtain uniprot id
  uniprot_id <- transcriptID2uniprotID(transcriptID)
  
  # obtain transcript length
  length <- transcriptID2length(transcriptID)
  
  # extract protien domain data
  protien_domain <- fetchDomain(uniprot_id)
  
  # construct gene from data collected
  geneData <- construct_gene(gene, protien_domain, length)
  
  # construct data frame of observed mutations
  observed_mutation <- mutationObs(data, fill_value, label_column, spread_degree)
  
  # construct data frame of cosmic mutations
  if(cosmic == TRUE)
  {
    cosmic_mutation <- mutationCOS(transcriptID, spread_degree)
  } else {
    cosmic_mutation <- NULL
  }
  
  # construct the lolliplot
  plot <- build_lolli(geneData, length, observed_mutation, cosmic_mutation, fill_value, label_column, plot_text_angle, plot_text_size, point_size, gene_colour)
  
  return(plot)
  
}