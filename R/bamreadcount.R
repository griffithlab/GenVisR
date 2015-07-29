#' Count nucleotide reads at SNP locations
#' 
#' Given the bam file path, count the number of reads at specified snp locations
#' @name bamreadcount
#' @param bamfile Path to the bam file
#' @param genome Object of class BSgenome corresponding to a genome of interest
#' @param targetbed Object of class data frame containing target locations in .bed format and containing columns chr, start, end
#' @return object of class data frame containing readcount information
#' @export 
#' @import GenomicRanges
#' @import Rsamtools
#' @importFrom "IRanges" IRanges
#' @import reshape2
#' @import Biostrings

bamreadcount <- function(bamfile, genome, targetbed = NULL)
{   
    # obtain the bam index file
    bai <- paste0(bamfile,".bai")
    
    # Perform basic quality checks on input data
    list <- bamreadcount.qual(bai, genome, targetbed)
    bai <- list[[1]]
    genome <- list[[2]]
    targetbed <- list[[3]]
    
    ## Read in the target bed locations. If none specified, read in 24 Pengelly loci.
    if(!is.null(targetbed)){
        pengelly <- read.table(file=targetbed, sep='\t', header = TRUE)
    }
    else{pengelly <- SNPloci}
    pengelly.chr <- pengelly
    
    ## Check if pengelly bed file has 'chr' associated with chromosome number in chr column.
    if(any(grepl("chr", pengelly$chr))){
        pengelly$chr <- gsub("chr","",pengelly$chr)
    }
    else{
        pengelly.chr$chr <- paste0("chr",pengelly.chr$chr)
    }
    grange <- with(pengelly, GRanges(chr, IRanges(start, end), strand=c('+')))
    grange.chr <- with(pengelly.chr, GRanges(chr, IRanges(start, end), strand=c('+')))
    
    ## Check using ScanBamHeader to see if 'chr' is included in chromosome names. Read in bam file
    what <- c("rname", "qname", "strand", "pos", "qwidth", "seq")
    chrcheck <- names(scanBamHeader(bamfile)[[1]]$targets[1])
    if(any(grepl("chr", chrcheck))){
        param <- ScanBamParam(which=grange.chr, what=what)
    }
    else{
        param <- ScanBamParam(which=grange, what=what)
    }
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
    x2 <- dcast(x, chr + position ~ nucleotide, fun.aggregate = sum)
    
    ## Add ref column and total_reads column
    getref <- as.character(getSeq(genome, grange.chr))
    total_reads <- data.frame(matrix(ncol=1,nrow=24))
    colnames(total_reads) <- 'total_reads'
    x3 <- cbind(x2[,1:2],getref,total_reads,x2[,3:6])
    ## total_reads column = sum of A, C, G, and T counts
    x3$total_reads <- x3$A+x3$C+x3$G+x3$T
    return(x3)
}