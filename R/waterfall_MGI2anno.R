#' Convert MGI File
#'
#' Convert columns of a mutation annotation file "MGI" into a format
#' recognizable by internal functions
#' @name waterfall_MGI2anno
#' @param x a data frame in MGI internal format
#' @param label_col Character string specifying the column name of a label
#' column
#' @return a data frame coerced from MGI to internal annotation format

waterfall_MGI2anno <- function(x, label_col, variant_class_order)
{
    # Check that correct column names are present and convert to internal format
    expec_col <- c('sample', 'gene_name', 'trv_type', 'chromosome_name',
                   'start', 'stop', 'reference', 'variant')
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
    x$key <- paste0(x$chromosome_name, ":", x$start, "-", x$stop, ":",
                    x$reference, "/", x$variant)
    x <- x[,c('key', 'sample', 'gene_name', 'trv_type', label_col)]
    
    if(!is.null(label_col))
    {
        colnames(x) <- c('key', 'sample', 'gene', 'trv_type', 'label')
    } else {
        colnames(x) <- c('key', 'sample', 'gene', 'trv_type')
    }
    
    # remove duplicate keys keeping entries based on a trv_type hiearchy
    if(is.null(variant_class_order)){
        mutation_order <- c("nonsense", "frame_shift_del",
                            "frame_shift_ins", "splice_site_del",
                            "splice_site_ins", "splice_site",
                            "nonstop", "in_frame_del", "in_frame_ins",
                            "missense", "splice_region_del",
                            "splice_region_ins", "splice_region",
                            "5_prime_flanking_region",
                            "3_prime_flanking_region",
                            "3_prime_untranslated_region",
                            "5_prime_untranslated_region", "rna",
                            "intronic", "silent", NA)
    } else {
        mutation_order <- variant_class_order
    }

    # Check that elements in trv_type are in the mutation order
    if(any(!x$trv_type %in% mutation_order))
    {
        memo <- paste0("Detected an invalid mutation type, valid values for ",
                       "MGI are: ", toString(mutation_order))
        stop(memo)
    }
    
    x$trv_type <- factor(x$trv_type, levels=mutation_order)
    x <- x[order(x$sample, x$key, x$gene),]
    x <- x[!duplicated(x[, c("sample", "key", "gene")]),]
    
    return(x)
}
