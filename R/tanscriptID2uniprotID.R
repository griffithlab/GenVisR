#' Convert trascript ID to uniprot ID
#' 
#' Convert an ensembl transcript ID to the corresponding uniprot ID
#' @name transcriptID2uniprotID
#' @param transcriptID String specifying ensembl transcript ID
#' @return String specifying uniprot ID

transcriptID2uniprotID <- function(transcriptID)
{
  ###############################################################################################
  ################### Function to retieve Uniprot ID for ensembl transcript #####################
  ###############################################################################################
  
  library("UniProt.ws")
  
  # select the database to grab the key from
  kt <- "ENSEMBL_TRANSCRIPT"
  
  # Query the transcript ID
  keys <- c(transcriptID)
  
  # Select data to retrieve back (uniprot ID)
  columns <- c("UNIPROTKB")
  
  # submit the query to return the uniprot ID
  result <- select(UniProt.ws, keys, columns, kt)
  uniprotID <- as.character(result[2])
  
  return(uniprotID)
  
}