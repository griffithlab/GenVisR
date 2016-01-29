#' Check input to TvTi
#'
#' Perform quality check for input to function TvTi
#' @name TvTi_qual
#' @param x Object of class data frame containing columns 'sample', reference',
#' 'variant' for 'MGI' file or 'Tumor_Sample_Barcode', 'Reference_Allele',
#' 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2' for 'MAF' file
#' @param y Object of class data frame containing columns "Prop", "trans_tranv"
#' @param z Object of class data frame containing columns "sample", "variable",
#' "value" denoting clinical information
#' @param file_type Character string spedifying th input file type expected
#' @return a data frame, or list of data frames passing quality checks

TvTi_qual <- function(x, y=NULL, z=NULL, file_type='MAF')
{
    # Check file type is valid
    if(!grepl("MAF|MGI", file_type))
    {
        memo <- paste0("Did not recognize input to paramter fileType as a",
                       " valid argument... Please specify one of \"MGI\"",
                       " or \"MAF\"")
        stop(memo)
    }
    
    # Check if x input is a data frame
    if(!is.data.frame(x))
    {
        memo <- paste0("argument supplied to x is not an object of class",
                       " data frame, attempting to coerce")
        warning(memo)
        x <- as.data.frame(x)
    }

    # check for duplicate elements in x
    if(nrow(unique(x)) != nrow(x))
    {
        warning("Detected duplicate rows in x, was this expected?")
    }

    # Check if y input is a data frame
    if(!is.null(y))
    {
        # Check y input if data frame
        if(is.data.frame(y))
        {
            if(!all(colnames(y) %in% c('Prop', 'trans_tranv')))
            {
                memo <- paste0("Did not detect correct column names in",
                               "input to y, missing one of \"Prop\",",
                               "\"trans_tranv\"")
                stop(memo)
            }
        }

        # Check y input if vector
        if(is.vector(y))
        {
            y <- as.data.frame(y)
            y$trans_tranv <- rownames(y)
            colnames(y) <- c('Prop', 'trans_tranv')

            if(typeof(y$Prop) != "double" & typeof(y$Prop) != "numeric")
            {
                stop("values found in y are not of type double or numeric")
            }
        }

        if(!is.data.frame(y))
        {
            memo <- paste0("input to y is not an object of class data frame",
                           " or named vector")
            stop(memo)
        }
    }

    # Check columns of x input and change to internal format
    if(file_type == 'MGI')
    {
        # Check that columns are named appropriatley, if not print error
        proper_names <- c("reference", "variant", "sample")
        if(all(proper_names %in% colnames(x)))
        {
            message("Found appropriate columns")
        } else {
            memo <- paste0("Could not find all columns requested, missing ",
                           "one of \"reference\", \"variant\", \"sample\"")
            stop(memo)
        }

        x <- x[,c('reference', 'variant', 'sample')]

    } else if(file_type == 'MAF') {
        proper_names <- c("Tumor_Sample_Barcode", "Reference_Allele",
                          "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
        if(all(proper_names %in% colnames(x)))
        {
            message("Found appropriate columns")
        } else {
            memo <- paste0("Could not find all columns requested, missing one ",
                           "of \"Tumor_Sample_Barcode\", \"Reference_Allele\",",
                           " \"Tumor_Seq_Allele1\", \"Tumor_Seq_Allele2\"")
            stop(memo)
        }
        # Convert MAF file to internal format
        x <- TvTi_convMaf(x)
    } else {
        memo <- paste0("TvTi requires a fileType specification, please",
                       "specify one of \"MAF\" or \"MGI\" based on the",
                       "argument supplied to parameter x. See docs for help.")
        stop(memo)
    }

    # Remove any indels present in the data
    x <- TvTi_rmIndel(x)
    # Warn about multi nucleotide codes
    x <- TvTi_rmMnuc(x)

    # Check that reference and variant columns only contain the proper codes
    ref_codes <- c('A', 'C', 'G', 'T', '-', 0)
    if(!all(toupper(x$reference) %in% toupper(ref_codes)))
    {
        memo <- paste0("Unrecognized Base Detected in reference column, ",
                       "expected values are: ", toString(ref_codes))
        stop(memo)
    } else if(!all(toupper(x$variant) %in% toupper(ref_codes))) {
        memo <- paste0("Unrecognized Base Detected in reference column, ",
                       "expected values are: ", toString(ref_codes))
        stop(memo)
    }

    # check y input for proper row names
    if(!is.null(y))
    {
        trans_tranv_names <- c("A->C or T->G (TV)", "A->G or T->C (TI)",
                               "A->T or T->A (TV)", "G->A or C->T (TI)",
                               "G->C or C->G (TV)", "G->T or C->A (TV)")
        if(!all(rownames(y) %in% trans_tranv_names))
        {
            memo <- paste0("Did not detect a value for all combinations of ",
                           "transitions/transversions, please specify input ",
                           "for each of the following levels: ",
                           toString(trans_tranv_names))
            stop(memo)
        }

        # check that y sums to 1 (i.e. its a true proportion among all elements)
        if(sum(y$Prop) != 1)
        {
            stop("The sum of elements in y should equal 1")
        }
    }
    
    # Check input data to clinDat
    if(!is.null(z))
    {
        if(!is.data.frame(z))
        {
            stop("Did not detect a data frame for input to clinDat")
        }
        z <- droplevels(z)
        
        if(!all(c('sample', 'variable', 'value') %in% colnames(z)))
        {
            stop("Did not detect correct sample names in clinDat")
        }
        
        if(!all(levels(x$sample) %in% levels(z$sample)))
        {
            memo <- paste0("Found a sample supplied to clinData not found",
                           " in the data frame supplied to x")
            warning(memo)
        }
    }

    return(list('input1'=x, 'input2'=y, 'input3'=z))
}
