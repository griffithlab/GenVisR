#' Count nucleotide reads at SNP locations
#' 
#' Given the bam file path, count the number of reads at the 24 SNP locations
#' @name bamreadcount.qual
#' @param bai Path to the bam index file
#' @return data frame
#' @import GenomicRanges
#' @import Rsamtools
#' @importFrom "IRanges" IRanges
#' @import reshape2
#' @import Biostrings

bamreadcount.qual<-function(bai, genome, targetbed){
    # Check to see if index file is found
    if(!file.exists(bai))
    {
        stop("Could not find bam index file")
    }
    
    # Check to see if genome is of class BSgenome
    if(!class(genome)[1]=="BSgenome")
    {
        stop("Genome must be of class BSgenome")
    }
    
    # 
    
    return(targetbed)
}