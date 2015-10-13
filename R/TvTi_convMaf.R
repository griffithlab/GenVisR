#' Convert .maf format to internal format
#'
#' Convert data frame in .maf format to an internally recogized format
#' @name TvTi_convMAF
#' @param x Object of class data frame containing columns
#' 'Tumor_Sample_Barcode', 'Reference_Allele' 'Tumor_Seq_Allele1',
#' 'Tumor_Seq_Allele2'
#' @return a data frame, with column names 'sample', 'reference', 'variant'

TvTi_convMaf <- function(x)
{
    # Take out the appropriate columns and format for each allele
    x <- x[,c('Tumor_Sample_Barcode', 'Reference_Allele', 'Tumor_Seq_Allele1',
    'Tumor_Seq_Allele2')]

    allele1 <- x[,c('Tumor_Sample_Barcode', 'Reference_Allele',
    'Tumor_Seq_Allele1')]
    colnames(allele1) <- c('sample', 'reference', 'variant')
    allele2 <- x[,c('Tumor_Sample_Barcode', 'Reference_Allele',
    'Tumor_Seq_Allele2')]
    colnames(allele2) <- c('sample', 'reference', 'variant')

    #!!!Developer note: The if's are here because subsetting when there is
    #!!! nothing to subset (integer0) causes problems

    # if the tumor allele 1 matchest tumor allele2 remove that information from
    # one of the alleles
    if(any(as.character(allele1$variant) == as.character(allele2$variant)))
    {
        allele1 <- allele1[-which(as.character(allele1$variant) == as.character(allele2$variant)),]
    }

    # if the allele matches the reference remove it from the data
    if(any(as.character(allele1$reference) == as.character(allele1$variant)))
    {
        allele1 <- allele1[-which(as.character(allele1$reference) == as.character(allele1$variant)),]
    }
    if(any(as.character(allele2$reference) == as.character(allele2$variant)))
    {
        allele2 <- allele2[-which(as.character(allele2$reference) == as.character(allele2$variant)),]
    }

    # bind the data from both alleles together
    x <- rbind(allele1, allele2)

    return(x)
}
