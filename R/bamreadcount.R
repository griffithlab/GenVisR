#' Count nucleotide reads at SNP locations
#' 
#' Given the bam file path, count the number of reads at the 24 SNP locations
#' @name bamreadcount
#' @param bamfile Path to the bam file
#' @return data frame
#' @export 
#' @import GenomicRanges
#' @import Rsamtools
#' @importFrom "IRanges" IRanges
#' @import reshape2

bamreadcount <- function(bamfile){
  ## Need index file
  bai <- paste(bamfile,".bai", sep='')
  
  ## Read in the 24 Pengelly SNP locations
  pengelly <- data.frame('chr'=c(1,1,2,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22))
  pengelly$start=pengelly$end <- c(67861520,179520506,169789016,227896976,4403767,5749904,82834630,146755140,48450157,94935937,100190780,100219314,16133413,993930,39433606,50769717,34528948,70303580,71197748,21413869,10267077,6100088,44323590,21141300)
  grange <- with(pengelly, GRanges(chr, IRanges(start, end)))
  
  ## Read in bam file (need index file name as well)
  what <- c("rname", "qname", "strand", "pos", "qwidth", "seq")
  param <- ScanBamParam(which=grange, what=what)
  bam <- scanBam(bamfile, bai, param=param)
  
  ## Pileup generates a table of nucleotide counts at each location by strand
  pileup_table <- pileup(bamfile, bai, scanBamParam = param)
  ## Remove strand and which_label columns
  pileup_table <- pileup_table[,c('seqnames','pos','nucleotide','count')]
  ## Rename columns to match bamreadcount table
  colnames(pileup_table)[1:2] <- c('chr','position')
  
  ## Use reshape to create separate columns for nucleotides
  ## Result (x2) is: chr, position, A, C, G, T
  x <- melt(pileup_table, c('chr','position','nucleotide'), 'count')
  x2 <- cast(x, chr + position ~ nucleotide, fun.aggregate = sum)
  
  ## Add ref column and total_reads column
  ref <- c('C','G','T','C','A','T','T','G','T','T','A','G','A','C','A','G','G','G','G','T','T','A','T','T')
  total_reads <- data.frame(matrix(ncol=1,nrow=24))
  colnames(total_reads) <- 'total_reads'
  x3 <- cbind(x2[,1:2],ref,total_reads,x2[,3:6])
  ## total_reads column = sum of A, C, G, and T counts
  x3$total_reads <- x3$A+x3$C+x3$G+x3$T
  return(x3)
}