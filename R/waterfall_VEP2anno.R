#' Convert VEP File
#'
#' Convert columns of a variant effect predictor file "VEP" into a format
#' recognizable by internal functions
#' @name waterfall_VEP2anno
#' @param x a data frame in VEP format
#' @param label_col Character string specifying the column name of a
#' label column
#' @return a data frame coerced from VEP to intermal expected format
#' @importFrom plyr rbind.fill

waterfall_VEP2anno <- function(x, label_col, variant_class_order)
{
    warning("VEP support is experimental, please double-check output!")
    
    # Make sure theere is an "Extra" Column
    if(!any(colnames(x) %in% c("Extra"))){
        memo <- paste("Did not detect a column name called \"Extra\",",
                      "the fields \"SYMBOL\" and \"IND\" within this",
                      "column is required!")
        stop(memo)
    }
    
    # Split fields in the "Extra" column of a VEP file into actual columns
    # and recombine with the original data in x
    a1 <- function(string){
        split_str <- do.call(rbind, strsplit(string, ";|="))
        name_indx <- seq(1, ncol(split_str), 2)
        value_indx <- seq(2, ncol(split_str), 2)
        str_dataframe <- as.data.frame(t(split_str[,value_indx]))
        colnames(str_dataframe) <- split_str[,name_indx]
        return(str_dataframe)
    }
    extra_col <- lapply(as.character(x$Extra), a1)
    extra_col <- do.call(plyr::rbind.fill, extra_col)
    x <- x[,-which(colnames(x) == "Extra")]
    x <- cbind(x, extra_col)
    
    # Check that correct column names are present and convert to internal format
    expec_col <- c('IND', 'SYMBOL', 'Consequence', 'Location', 'Allele')
    
    if(!is.null(label_col))
    {
        expec_col <- c(expec_col, label_col)
    }
    
    if(!all(expec_col %in% colnames(x)))
    {
        memo <- paste0("Did not detect correct column names, column names",
                       " should be: ", toString(expec_col), " these should be",
                       " either column names or fields in the Extra column!")
        stop(memo)
    }
    
    # add a unique key and remove uneccessary columns
    x$key <- paste0(x$Location, ":", x$Allele)
    x <- x[,c('key', 'IND', 'SYMBOL', 'Consequence', label_col)]
    
    if(!is.null(label_col))
    {
        colnames(x) <- c('key', 'sample', 'gene', 'trv_type', 'label')
    } else {
        colnames(x) <- c('key', 'sample', 'gene', 'trv_type')
    }
    
    # Take the first element in the consequence field (now trv_type) and remove
    # the rest (1st element should have highest impact in vep)
    a2 <- function(string){
        string <- as.character(string)
        splitStr <- strsplit(string, ",")
        return(splitStr[[1]][1])
    }
    x$trv_type <- sapply(x$trv_type, a2)
    
    # Remove any NA values in the data
    if(any(is.na(x$gene))){
        gene_row_count <- nrow(x)
        x <- x[which(is.na(x$gene)),]
        
        memo <- paste("Removed", gene_row_count - nrow(x), "entries without a",
                      "gene symbol!")
        message(memo)
    }
    
    if(any(is.na(x$sample))){
        sample_row_count <- nrow(x)
        x <- x[which(is.na(x$sample)),]
        
        memo <- paste("Removed", sample_row_count - nrow(x), "entries without",
                      "an annotated sample!")
        message(memo)
    }
    
    # remove duplicate keys keeping entries based on a trv_type hiearchy
    if(is.null(variant_class_order)){
        mutation_order <- c("transcript_ablation", "splice_acceptor_variant",
                            "splice_donor_variant", "stop_gained",
                            "frameshift_variant", "stop_lost", "start_lost",
                            "transcript_amplification", "inframe_insertion",
                            "inframe_deletion", "missense_variant",
                            "protein_altering_variant", 
                            "splice_region_variant",
                            "incomplete_terminal_codon_variant",
                            "stop_retained_variant", "synonymous_variant",
                            "coding_sequence_variant",
                            "mature_miRNA_variant", "5_prime_UTR_variant",
                            "3_prime_UTR_variant",
                            "non_coding_transcript_exon_variant",
                            "intron_variant", "NMD_transcript_variant",
                            "non_coding_transcript_variant",
                            "upstream_gene_variant",
                            "downstream_gene_variant", "TFBS_ablation",
                            "TFBS_amplification", "TF_binding_site_variant",
                            "regulatory_region_ablation",
                            "regulatory_region_amplification",
                            "feature_elongation",
                            "regulatory_region_variant",
                            "feature_truncation", "intergenic_variant")
    } else {
        mutation_order <- variant_class_order
    }
    
    # Check that elements in trv_type are in the mutation order
    if(any(!x$trv_type %in% mutation_order))
    {
        memo <- paste0("Detected an invalid mutation type, valid values for ",
                       "VEP are: ", toString(mutation_order))
        stop(memo)
    }

    x$trv_type <- factor(x$trv_type, levels=mutation_order)
    x <- x[order(x$sample, x$key, x$gene),]
    x <- x[!duplicated(x[, c("sample", "key", "gene")]),]
    
    return(x)
}