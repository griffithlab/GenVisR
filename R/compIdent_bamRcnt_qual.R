#' Count nucleotide reads at SNP locations
#'
#' Given the bam file path, count the number of reads at the 24 SNP locations
#' @name compIdent_bamRcnt_qual
#' @param genome Object of class BSgenome corresponding to a genome of interest
#' @param targetbed Object of class data frame containing target locations in .bed format
#' @return list of data objects passing quality checks
#' @noRd

compIdent_bamRcnt_qual<-function(genome, targetbed)
{
    # Check to see if genome is of class BSgenome
    if(!is(genome, "BSgenome"))
    {
        stop("Genome must be of class BSgenome")
    }

    # Check bed file input
    if(!is.null(targetbed))
    {
        if(!is.data.frame(targetbed))
        {
            memo <- paste0("Object supplied to target bed does not appear to",
                           " be a data frame... attempting to coerce")
            warning(memo)
            targetbed <- as.data.frame(targetbed)
            colnames(targetbed) <- c("chr", "start", "end")
        }

        if(!all(c('chr', 'start', 'end') %in% colnames(targetbed)))
        {
            memo <- paste0("Did not detect correct columns in targetbed, ",
                           "missing one of \"chr\", \"start\", \"end\"")
            stop(memo)
        }
    }

    return(list(genome, targetbed))
}
