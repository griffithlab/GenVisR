#' Count nucleotide reads at SNP locations
#' 
#' Given the bam file path, count the number of reads at the 24 SNP locations
#' @name bamreadcount.qual
#' @param bai Path to the bam index file
#' @param genome Object of class BSgenome corresponding to a genome of interest
#' @param targetbed Object of class data frame containing target locations in .bed format
#' @return list of data objects passing quality checks

bamreadcount.qual<-function(bai, genome, targetbed)
{
    # Check to see if index file is found
    if(!file.exists(bai))
    {
        stop("Could not find bam index file for: ", gsub(".bai$", "bam", bai))
    }
    
    # Check to see if genome is of class BSgenome
    if(!class(genome)[1]=="BSgenome")
    {
        stop("Genome must be of class BSgenome")
    }
    
    # Check bed file input
    if(!is.null(targetbed))
    {
        if(!is.data.frame(targetbed))
        {
            warning("Object supplied to target bed does not appear to be a data frame, attempting to coerce...")
            targetbed <- as.data.frame(targetbed)
            colnames(targetbed) <- c("chr", "start", "end")
        }
        
        if(!all(c('chr', 'start', 'end') %in% colnames(targetbed)))
        {
            stop("Did not detect correct columns in targetbed, missing one of chr, start, end")
        }
    }
    
    return(list(bai, genome, targetbed))
}