#' Convert .maf format to internal format
#' 
#' Convert data frame in .maf format to an internally recogized format
#' @name TvTi.convMAF
#' @param x Object of class data frame containing columns 'Tumor_Sample_Barcode', 'Reference_Allele' 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2'
#' @return a data frame, with column names 'sample', 'reference', 'variant'

TvTi.convMaf <- function(x)
{
  # Take out the appropriate columns and format for each allele
  x <- x[,c('Tumor_Sample_Barcode', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2')]
  
  allele1 <- x[,c('Tumor_Sample_Barcode', 'Reference_Allele', 'Tumor_Seq_Allele1')]
  colnames(allele1) <- c('sample', 'reference', 'variant')
  allele2 <- x[,c('Tumor_Sample_Barcode', 'Reference_Allele', 'Tumor_Seq_Allele2')]
  colnames(allele2) <- c('sample', 'reference', 'variant')
  
  # bind the data together
  x <- unique(rbind(allele1, allele2))
  
  return(x)
}