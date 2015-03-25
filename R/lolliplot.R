#' Construct Lolliplot from data
#' 
#' Construct a Lolliplot from object of class data frame giving observed mutations and transcript
#' @name lolliplot
#' @param data object of class data frame containing columns transcript_name, gene, and amino_acid_change for specific transcript
#' @param cosmic boolean value specifying if cosmic mutations should be retrieved and plotted
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