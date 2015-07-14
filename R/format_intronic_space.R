#' format intronic space
#' 
#' format intronic space
#' @name format_intronic_space
#' @param introns list of data frames containing intronic features
#' @param exons list of data frames contatining exonic/UTR features
#' @param gr Granges object specifying the ROI
#' @return Object of class data frame

format_intronic_space <- function(introns, exons, gr)
{
  # Determine what chromosome the GRange object refers to
  chr <- unique(as.character(seqnames(gr)))
  
  # bind intronic and exonic data together
  gene <- lapply(exons, rbind, introns)
  
  # Calculate any pseudo_exons present (i.e. a genomic coordinate present as an exon/UTR in one isoform but not another)
  gene <- lapply(gene, Calculate_Pseudoexons, chr)
  
  # compute the transform
  gene <- lapply(gene, intronic_transform)
  
  return(gene)
}