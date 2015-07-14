#' Convert trascript ID to protien length
#' 
#' Retrieve the length of an ensembl transcript ID in amino acids
#' @name transcriptID2length
#' @param transcriptID String specifying ensembl transcript ID
#' @param up uniprot data object from UniProt.ws
#' @return numeric value specifying protien length
#' @import UniProt.ws

transcriptID2length <- function(transcriptID, up)
{
  ##############################################################################################
  #################### Function to retrieve ensembl transcript amino acid length ###############
  ##############################################################################################

  # select the database to grab the key from
  kt <- "ENSEMBL_TRANSCRIPT"
  
  # query the transcript ID
  keys <- c(transcriptID)
  
  # select data to retrieve back (AA length)
  columns <- c("LENGTH")
  
  # submit the query to return the length in amino acid residues
  result <- select(up, keys, columns, kt)
  AAlength <- as.numeric(result[2])
  
  return(AAlength)
}