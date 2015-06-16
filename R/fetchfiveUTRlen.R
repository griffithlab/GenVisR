#' fetch 5'UTR length
#' 
#' Retrieve 5'UTR length for a given ensembl transcript
#' @name fetchfiveUTRlen
#' @param enstID a character string specifying the ensembl transcript ID
#' @param dataset a character string specifying the name of the biomart dataset to use
#' @return an integer corresponding to the 5'UTR length
#' @import biomaRt

fetchfiveUTRlen <- function(enstID, dataset="hsapiens_gene_ensembl")
{
  ##############################################################################################
  #################### Function to retrieve protien domains given uniprot id ###################
  ##############################################################################################
  
  # Load in database and select dataset
  ensembl_mart <- useMart("ensembl")
  ensembl_mart <- useDataset(dataset, mart=ensembl_mart)
  
  # Print message
  datalist <- listDatasets(ensembl_mart)
  assembly <- datalist[which(datalist$dataset == dataset),c("description")]
  message("using ", assembly, " to convert notation")
  warning("output may be incorrect if assembly for amino acid change differs from ", assembly)
  warning("c.notation conversion is experimental, this is not advised, check output!!!!!!!!")
  
  # Select attributes to retrieve 5'UTR (start, stop)
  attributes <- c("5_utr_start", "5_utr_end")
  
  # Apply various filters using vector of values
  filters <- c("ensembl_transcript_id")
  values <- c(as.character(enstID))
  
  # Retrieve data and calculate the total 5'UTR length
  FiveUTR <- getBM(attributes=attributes, filters=filters, values=values, mart=ensembl_mart)
  FiveUTR <- sum(FiveUTR[,2] - FiveUTR[,1], na.rm=TRUE)
  
  return(FiveUTR)
}