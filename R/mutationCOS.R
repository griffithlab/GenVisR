#' retrieve cosmic mutation observations
#' 
#' Create a data frame of cosmic mutation observations
#' @name mutationCOS
#' @param transcriptID character string specifying ensembl transcript ID to retrieve mutations for
#' @return object of class data frame giving comsic mutation observations
#' @import biomaRt

mutationCOS <- function(transcriptID, rep.fact, rep.dist.lmt, attr.fact, adj.max, adj.lmt, iter.max)
{
  #################################################################################
  ################## function to return cosmic variants given transcript id #######
  #################################################################################

  # Load in database and select dataset
  ensembl_mart <- useMart("ensembl")
  ensembl_mart <- useDataset("hsapiens_gene_ensembl", mart=ensembl_mart)
  
  # Select attributes to retrieve (protein domain, start, stop)
  attributes <- c("somatic_variation_name", "somatic_peptide_location")
  
  # Apply various filters using vector of values
  filters <- c("ens_hs_transcript", "somatic_variation_source", "with_validated_snp")
  values <- list(transcriptID, "COSMIC", TRUE)
  
  # Retrieve data
  message("Retrieving cosmic mutations")
  cosmic <- getBM(attributes=attributes, filters=filters, values=values, mart=ensembl_mart)
  if(nrow(cosmic) == 0)
  {
    stop("No Cosmic Mutations retrieved, check biomart query in mutationCos.R")
  }
  
  # omit values with NA at position (need to look at why this is ocurring)
  cosmic <- na.omit(cosmic)
  
  # re-label columns and add another column for y axis
  colnames(cosmic) <- c("cosmic_id", "mutation_coord")
  cosmic$height_min <- -1
  
  # Dodge mutation coordinates on the x axis
  message("applying force field to cosmic mutations")
  cosmic <- cosmic[order(cosmic$mutation_coord),] 
  cosmic$coord_x_dodge <- dodge_coord_x(as.vector(cosmic$mutation_coord), rep.fact=rep.fact, rep.dist.lmt=rep.dist.lmt, attr.fact=attr.fact, adj.max=adj.max, adj.lmt=adj.lmt, iter.max=iter.max)
  
  # Redefine and return grouping information and then dodge y coordinates
  cosmic$group <- group_mutation_coord(as.vector(cosmic$mutation_coord))
  cosmic$coord_y_dodge <- dodge_coord_y(cosmic, track='bottom')
  
  return(cosmic)
}