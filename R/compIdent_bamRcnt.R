#' Count nucleotide reads at SNP locations
#'
#' Given the bam file path, count the number of reads at specified snp locations
#' @name compIdent_bamRcnt
#' @param bamfile Path to the bam file
#' @param genome Object of class BSgenome corresponding to a genome of interest
#' @param target Object of class data frame containing target locations in
#' 1-base format and containing columns "chr", "start", "end", "var", "name"
#' @param debug Boolean specifying if test datasets should be used for
#' debugging.
#' @return object of class data frame containing readcount information
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBam
#' @importFrom Rsamtools pileup
#' @importFrom data.table melt
#' @importFrom data.table dcast
#' @noRd

compIdent_bamRcnt <- function(bamfile, genome, target=NULL, debug=FALSE)
{

    # Perform basic quality checks on input data
    list <- compIdent_bamRcnt_qual(genome, target)
    genome <- list[[1]]
    target <- list[[2]]

    # Read in the target bed locations. If none specified, read in 24 Pengelly
    # loci.
    if(!is.null(target))
    {
        target <- target
    } else {
        target <- GenVisR::SNPloci
    }
    
    # Target locations must conform to the bam files being read in, create
    # two versions one with "chr1" and the other with "1" and use whichever is
    # appropriate
    if(any(grepl("chr", target$chr)))
    {
        target$chr <- gsub("chr", "", target$chr)
    } else {
        target.chr <- target
        target.chr$chr <- paste0("chr", target$chr)
    }
    
    # Similar to above create two Grange objects one with "chr" appendix and
    # one without
    grange <- GenomicRanges::GRanges(target$chr,
                                     IRanges::IRanges(target$start,
                                                      target$end))
    
    grange.chr <- GenomicRanges::GRanges(target.chr$chr,
                                         IRanges::IRanges(target.chr$start,
                                                          target.chr$end))
    
    # If debug flag is true run the test case else read in specified bam
    # files
    if(isTRUE(debug))
    {
        pileup_table <- bamfile
    } else{
        # obtain the bam index file name
        bai <- paste0(bamfile,".bai")
        
        # Check to see if index file is found
        if(!file.exists(bai))
        {
            memo <- paste0("Could not find bam index file for: ",
                           gsub(".bai$", "bam", bai))
            stop(memo)
        }
        
        # Check using ScanBamHeader to see if 'chr' is included in chromosome
        # names
        what <- c("rname", "qname", "pos", "qwidth", "seq")
        chrcheck <- names(Rsamtools::scanBamHeader(bamfile)[[1]]$targets[1])
        
        # set up the appropriate param based on if chr is present in bam or not
        if(any(grepl("chr", chrcheck)))
        {
            param <- Rsamtools::ScanBamParam(which=grange.chr, what=what)
        } else {
            param <- Rsamtools::ScanBamParam(which=grange, what=what)
        }
        
        # Pileup generates a table of nucleotide counts at each location
        pileup_table <- Rsamtools::pileup(bamfile, bai, scanBamParam=param)
    }
    
    # Remove strand and which_label columns
    pileup_table <- pileup_table[, !names(pileup_table) %in% c('strand',
                                                               'which_label')]
    
    # Rename columns to match bamreadcount table
    colnames(pileup_table)[1:2] <- c('chr','position')

    # Use reshape to create separate columns for nucleotides
    # Result (x2) is: chr, position, A, C, G, T
    x <- data.table::melt(pileup_table, c('chr','position','nucleotide'), 'count')
    x <- data.table::dcast(x, chr + position ~ nucleotide, fun.aggregate = sum)

    # Add ref column
    x$getref <- as.character(Biostrings::getSeq(genome, grange.chr))
    
    # Add total reads columns
    x$total_reads <- rowSums(x[,c('A', 'C', 'G', 'T')])
    
    # Add back in var and name colummns
    x$key <- paste0(x$chr, ":", x$position)
    target$key <- paste0(target$chr, ":", target$end)
    x <- merge(x, target, by=c("key", "chr"))
    x <- x[,c("chr", "position", "A", "C", "G", "T", "total_reads", "getref",
              "var", "name")]
    
    return(x)
}