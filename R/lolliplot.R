#' Construct Lolliplot from data
#' 
#' Construct a Lolliplot from object of class data frame giving observed mutations and an ensembl transcript id
#' @name lolliplot
#' @param x object of class data frame containing columns transcript_name, gene, and amino_acid_change and rows denoting mutations
#' @param cosmic boolean value specifying if cosmic mutations should be retrieved and plotted
#' @param fill_value character string giving the name of the column to shade variants on "required"
#' @param label_column character string specifying column containing text information to be plotted, defaults to NULL
#' @param plot_text_angle integer specifying angle of text to be plotted if label_column is specified
#' @param plot_text_size integer specifying the size of the text plotted if label_column is specified
#' @param point_size integer specifying the size of points plotted
#' @param gene_colour color specifying the background fill of the plotted gene
#' @param obs.rep.fact repulsive factor for plotted mutations observed track
#' @param obs.rep.dist.lmt repulsive distance limit for plotted mutations observed track
#' @param obs.attr.fact attraction factor for plotted mutations observed track
#' @param obs.adj.max maximum position change for each iteration observed track
#' @param obs.adj.lmt position adjustment limit which simulation stops observed track
#' @param obs.iter.max maximum iterations beyond which to stop the simulation observed track
#' @param cos.rep.fact repulsive factor for plotted mutations cosmic track
#' @param cos.rep.dist.lmt repulsive distance limit for plotted mutations cosmic track
#' @param cos.attr.fact attraction factor for plotted mutations cosmic track
#' @param cos.adj.max maximum position change for each iteration cosmic track
#' @param cos.adj.lmt position adjustment limit which simulation stops cosmic track
#' @param cos.iter.max maximum iterations beyond which to stop the simulation cosmic track
#' @param plot_sidechain boolean value whether to plot the amino acid sidechains instead of protein domains
#' @param taxId integer specifying the uniprot taxonomy id for the species of interest
#' @param ensembl.dataset character string specifying what ensembl data set to query when converting to p. notation, ignored if data is already in p.notation
#' @return object of class ggplot2
#' @export
#' @import UniProt.ws
#' @import RCurl

lolliplot <- function(x, cosmic=FALSE, fill_value=NULL, label_column=NULL, plot_text_angle=45, plot_text_size=5, point_size=3, gene_colour='#999999', obs.rep.fact=5000, obs.rep.dist.lmt=500, obs.attr.fact=.1, obs.adj.max=.1, obs.adj.lmt=.5, obs.iter.max=50000, cos.rep.fact=5000, cos.rep.dist.lmt=500, cos.attr.fact=.1, cos.adj.max=.1, cos.adj.lmt=.5, cos.iter.max=50000, plot_sidechain=FALSE, taxId=9606, ensembl.dataset="hsapiens_gene_ensembl")
{  
  # Perform quality check
  x <- lolliplot.qual(x)
  
  # Define a taxonomy ID for use in the "transcriptID2" function family for use with UniProt.ws
  up <- UniProt.ws(taxId=taxId)
  
  # extract transcript id
  transcriptID <- as.character(x$transcript_name[1])
  
  # extract HUGO gene name
  gene <- as.character(x$gene[1])
  
  # obtain uniprot id
  uniprot_id <- transcriptID2uniprotID(transcriptID, up)
  
  # obtain transcript length
  length <- transcriptID2length(transcriptID, up)
  
  # obtain amino acid sequence and format if it is requested to plot the sidechain
  if(plot_sidechain==T)
  {
    AAsequence <- transcriptID2sequence(transcriptID, up)
    AAsequence$sidechain <- sapply(AAsequence[,1], AA2sidechain)
  } else {
    AAsequence <- NULL
  }

  # extract protien domain data
  protien_domain <- fetchDomain(uniprot_id)
  
  # construct gene from data collected
  geneData <- construct_gene(gene, protien_domain, length)
  
  # construct data frame of observed mutations
  observed_mutation <- mutationObs(x, fill_value, label_column, obs.rep.fact, obs.rep.dist.lmt, obs.attr.fact, obs.adj.max, obs.adj.lmt, obs.iter.max, ensembl.dataset=ensembl.dataset)
  
  # construct data frame of cosmic mutations
  if(cosmic == TRUE)
  {
    cosmic_mutation <- mutationCOS(transcriptID, cos.rep.fact, cos.rep.dist.lmt, cos.attr.fact, cos.adj.max, cos.adj.lmt, cos.iter.max)
  } else {
    cosmic_mutation <- NULL
  }
  
  # construct the lolliplot
  plot <- build_lolli(geneData, length, observed_mutation, cosmic_mutation, fill_value, label_column, plot_text_angle, plot_text_size, point_size, gene_colour, AAsequence, plot_sidechain=plot_sidechain)
  
  return(plot)
}