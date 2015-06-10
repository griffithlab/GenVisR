#' fetch protein domains
#' 
#' Retrieve protein domains from interprot database given uniprot ID
#' @name fetchDomain
#' @param uniprotID String specifying uniprot ID
#' @return data frame of protien domains and start/stop coordinates
#' @import biomaRt

fetchDomain <- function(uniprotID)
{
  ##############################################################################################
  #################### Function to retrieve protien domains given uniprot id ###################
  ##############################################################################################

  # Load in database and select dataset
  interpro_mart <- useMart("prod-intermart_1")
  interpro_mart <- useDataset("protein", mart=interpro_mart)
  
  # Select attributes to retrieve (protein domain, start, stop)
  attributes <- c("entry_name", "pos_from", "pos_to")
  
  # Apply various filters using vector of values
  filters <- c("entry_type", "protein_accession", "fragment", "method_database_name")
  values <- as.list(c("Domain", uniprotID, "N", "Pfam"))
  
  # Retrieve data
  domains <- getBM(attributes=attributes, filters=filters, values=values, mart=interpro_mart)
  
  return(domains)
}