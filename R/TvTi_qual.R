#' Check input to TvTi
#' 
#' Perform quality check for input to function TvTi
#' @name TvTi_qual
#' @param x Object of class data frame containing columns 'sample', reference',
#' 'variant' for 'MGI' file or 'Tumor_Sample_Barcode', 'Reference_Allele',
#' 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2' for 'MAF' file
#' @param y Object of class data frame containing columns "Prop", "trans_tranv"
#' @param file_type Character string spedifying th input file type expected
#' @return a data frame, or list of data frames passing quality checks

TvTi_qual <- function(x, y=NULL, file_type='MAF')
{
    # Check if x input is a data frame
    if(!is.data.frame(x))
    {
        warning(x, "is not an object of class data frame, attempting to coerce")
        x <- as.data.frame(x)
        x <- relevel(x)
    }
    
    # check for duplicate elements in x
    if(nrow(unique(x)) != nrow(x))
    {
        warning("Detected duplicate rows in x, was this expected?")
    }
    
    # Check if y input is a data frame 
    if(!is.null(y))
    {
        # Check y input
        if(is.data.frame(y))
        {
            if(colnames(y) %in% c('Prop', 'trans_tranv'))
            {
                memo <- paste0("Did not detect correct columns names in",
                               "input to y, missing one of \"Prop\",",
                               "\"trans_tranv\"")
                stop(memo)
            }
        }
        
        if(is.vector(y))
        {
            y <- as.data.frame(y)
            y$trans_tranv <- rownames(y)
            colnames(y) <- c('Prop', 'trans_tranv')
            
            if(typeof(y$Prop) != "double" & typeof(y$Prop) != "numeric")
            {
                stop("values in", y, "are not of type double or numeric")
            }
        }
        
        if(!is.data.frame(y))
        {
            memo <- paste0(y, " is not an object of class data frame",
                           " attempting to coerce")
            warning(memo)
            y <- as.data.frame(y)
        }
    }
    
    # Check columns of x input and change to internal format
    if(file_type == 'MGI')
    {
        # Check that columns are named appropriatley, if not print error
        if(any(grepl('^reference$', colnames(x))) &&
               any(grepl('^variant$', colnames(x))) &&
               any(grepl('^sample$', colnames(x))))
        {
            message("Found appropriate columns")
        } else {
            memo <- paste0("Could not find all columns requested, missing ", 
                           "one of \"reference\", \"variant\", \"sample\"")
            stop(memo)
        } 
        
        x <- x[,c('reference', 'variant', 'sample')]
        
    } else if(file_type == 'MAF') {
        if(any(grepl('^Tumor_Sample_Barcode$', colnames(x))) &&
               any(grepl('^Reference_Allele$', colnames(x))) &&
               any(grepl('^Tumor_Seq_Allele1$', colnames(x))) &&
               any(grepl('^Tumor_Seq_Allele2$', colnames(x))))
        {
            message("Found appropriate columns")
        } else {
            memo <- paste0("Could not find all columns requeste, missing one ",
                           "of \"Tumor_Sample_Barcode\", \"Reference_Allele\",",
                           " \"Tumor_Seq_Allele1\", \"Tumor_Seq_Allele2\"")
            stop(memo)
        }
        # Convert MAF file to internal format
        x <- TvTi_convMaf(x)
    }
    
    # Remove any indels present in the data
    x <- TvTi_rmIndel(x)
    # Warn about multi nucleotide codes
    x <- TvTi_rmMnuc(x)
    
    # Check that reference and variant columns only contain the proper codes
    ref_codes <- c('A', 'C', 'G', 'T', '-', 0)
    if(!all(toupper(x$reference) %in% toupper(ref_codes)))
    {
        stop("Unrecognized Base Detected in reference column")
    } else if(!all(toupper(x$variant) %in% toupper(ref_codes))) {
        stop("Unrecognized Base Detected in variant column")
    }
    
    # check y input for proper row names
    if(!is.null(y))  
    {
        trans_tranv_names <- c("A->C or T->G", "A->G or T->C", "A->T or T->A",
                               "G->A or C->T", "G->C or C->G", "G->T or C->A")
        if(!all(rownames(y) %in% trans_tranv_names))
        {
            stop("Did not detect correct names in y")
        }
        
        # check that y sums to 1 (i.e. its a true proportion among all elements)
        if(sum(y$Prop) != 1)
        {
            stop("The sum of elements in y should equal 1")
        }
    }
    
    return(list('input1'=x, 'input2'=y))
}
