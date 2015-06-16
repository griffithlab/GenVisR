#' Check input to lolliplot
#' 
#' Perform Basic quality checks for lolliplot input
#' @name lolliplot.qual
#' @param data object of class data frame containing columns transcript_name, gene, and amino_acid_change and rows denoting mutations
#' @return lolliplot data frame passing basic quality control

lolliplot.qual <- function(x)
{
  # Check for internet connectivity
  if(!is.character(getURL("www.google.com")))
  {
    stop("Did not detect an internet connection, check internet connectivity")
  }
  
  # Check for correct columns
  if(!all(c('transcript_name', 'gene', 'amino_acid_change') %in% colnames(x)))
  {
    Stop("Did not detect correct columns, missing one of transcript_name, gene, amino_acid_change")
  }
  
  return(x)
}