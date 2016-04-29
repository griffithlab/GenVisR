#' Convert MAF File
#'
#' Convert columns of a mutation annotation file "MAF" into a format
#' recognizable by internal functions
#' @name waterfall_MAF2anno
#' @param x a data frame in MAF format
#' @param label_col Character string specifying the column name of a
#' label column
#' @return a data frame coerced from MAF to TGI format

waterfall_MAF2anno <- function(x, label_col, variant_class_order)
{
    # Check that correct column names are present and convert to internal format
    expec_col <- c('Chromosome', 'Start_Position', 'End_Position',
                   'Reference_Allele', 'Tumor_Seq_Allele2',
                   'Tumor_Sample_Barcode', 'Hugo_Symbol',
                   'Variant_Classification')

    if(!is.null(label_col))
    {
        expec_col <- c(expec_col, label_col)
    }

    if(!all(expec_col %in% colnames(x)))
    {
        memo <- paste0("Did not detect correct column names, column names
                       should be: ", toString(expec_col))
        stop(memo)
    }

    # add a unique key and remove uneccessary columns
    x$key <- paste0(x$Chromosome, ":", x$Start_Position, "-", x$End_Position,
                    ":", x$Reference_Allele, "/", x$Tumor_Seq_Allele2)
    x <- x[,c('key', 'Tumor_Sample_Barcode', 'Hugo_Symbol',
              'Variant_Classification', label_col)]

    if(!is.null(label_col))
    {
        colnames(x) <- c('key', 'sample', 'gene', 'trv_type', 'label')
    } else {
        colnames(x) <- c('key', 'sample', 'gene', 'trv_type')
    }
    
    # remove duplicate keys keeping entries based on a trv_type hiearchy
    if(is.null(variant_class_order)){
        mutation_order <- c("Nonsense_Mutation", "Frame_Shift_Ins",
                            "Frame_Shift_Del", "Translation_Start_Site",
                            "Splice_Site", "Nonstop_Mutation",
                            "In_Frame_Ins", "In_Frame_Del",
                            "Missense_Mutation", "5\'Flank",
                            "3\'Flank", "5\'UTR", "3\'UTR", "RNA", "Intron",
                            "IGR", "Silent", "Targeted_Region", NA)
    } else {
        mutation_order <- variant_class_order
    }
    
    # Check that elements in trv_type are in the mutation order
    if(any(!x$trv_type %in% mutation_order))
    {
        memo <- paste0("Detected an invalid mutation type, valid values for ",
                       "MAF are: ", toString(mutation_order))
        stop(memo)
    }

    x$trv_type <- factor(x$trv_type, levels=mutation_order)
    x <- x[order(x$sample, x$key, x$gene),]
    x <- x[!duplicated(x[, c("sample", "key", "gene")]),]
    
    return(x)
}
