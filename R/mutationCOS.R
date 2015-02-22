#' retrieve cosmic mutation observations
#' 
#' Create a data frame of cosmic mutation observations
#' @name mutationCOS
#' @param transcriptID character string specifying ensembl transcript ID to retrieve mutations for
#' @return object of class data frame giving comsic mutation observations

mutationCOS <- function(transcriptID)
{
  #################################################################################
  ################## function to return cosmic variants given transcript id #######
  #################################################################################
  
  library("biomaRt")
  
  # Load in database and select dataset
  ensembl_mart <- useMart("ensembl")
  ensembl_mart <- useDataset("hsapiens_gene_ensembl", mart=ensembl_mart)
  
  # Select attributes to retrieve (protein domain, start, stop)
  attributes <- c("somatic_variation_name", "somatic_peptide_location")
  
  # Apply various filters using vector of values
  filters <- c("ens_hs_transcript", "somatic_variation_source", "with_validated_snp")
  values <- list(transcriptID, "COSMIC", TRUE)
  
  # Retrieve data
  cosmic <- getBM(attributes=attributes, filters=filters, values=values, mart=ensembl_mart)
  
  # omit values with NA at position (need to look at why this is ocurring)
  cosmic <- na.omit(cosmic)
  
  # re-label columns and add another column for y axis
  colnames(cosmic) <- c("cosmic_id", "mutation_coord")
  cosmic$height_min <- -1.5
  
  return(cosmic)
}