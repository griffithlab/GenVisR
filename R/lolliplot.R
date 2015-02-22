#' Construct Lolliplot from data
#' 
#' Construct a Lolliplot from object of class data frame giving observed mutations and transcript
#' @name lolliplot
#' @param data object of class data frame containing columns transcript_name, gene, and amino_acid_change for specific transcript
#' @param cosmic boolean value specifying if cosmic mutations should be retrieved and plotted
#' @return object of class ggplot2
#' @export

lolliplot <- function(data, cosmic=TRUE)
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
  observed_mutation <- mutationObs(data)
  
  # construct data frame of cosmic mutations
  if(cosmic == TRUE)
  {
    cosmic_mutation <- mutationCOS(transcriptID)
  } else {
    cosmic_mutation <- NULL
  }
  
  # construct the lolliplot
  plot <- build_lolli(geneData, length, observed_mutation, cosmic_mutation)
  
  return(plot)
  
}