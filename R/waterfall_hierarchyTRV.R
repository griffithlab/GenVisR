#' Hiearchical removal of MAF entries
#'
#' Remove MAF entries with the same gene/sample in an ordered fashion such that
#' the most deleterious are retained
#' @name waterfall_hierarchyTRV
#' @param x a data frame in long format with columns sample, gene,
#' trv_type
#' @param file_type The type of file to act on one of 'MAF", "MGI", "Custom"
#' @param variant_class_order character vector giving the hierarchical order of
#' mutation types to plot
#' @return a data frame with multiple mutations in the same sample/gene
#' collapsed on the most deleterious

waterfall_hierarchyTRV <- function(x, file_type, variant_class_order)
{
    message("setting mutation hierarchy...")
    # if variant_class_order is null use predefined values
    if(is.null(variant_class_order))
    {
        if(toupper(file_type) == toupper('MGI'))
        {
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
        } else if(toupper(file_type) == toupper('MAF')) {
            mutation_order <- c("Nonsense_Mutation", "Frame_Shift_Ins",
                                "Frame_Shift_Del", "Translation_Start_Site",
                                "Splice_Site", "Nonstop_Mutation",
                                "In_Frame_Ins", "In_Frame_Del",
                                "Missense_Mutation", "5\'Flank",
                                "3\'Flank", "5\'UTR", "3\'UTR", "RNA", "Intron",
                                "IGR", "Silent", "Targeted_Region", NA)
        } else if(toupper(file_type) == toupper('VEP')) {
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
        } else if(toupper(file_type) == toupper('Custom')) {
            memo <- paste0("Detected NULL in variant_class_order, ",
                           "this parameter is required if file_type is set ",
                           "to \"Custom\"")
            stop(memo)
        }
    } else {
        mutation_order <- variant_class_order
    }

    # Check that elements in trv_type are in the mutation order
    if(any(!x$trv_type %in% mutation_order))
    {
        memo <- paste0("Detected an invalid mutation type, valid values for ",
                       file_type, " are: ", toString(mutation_order))
        stop(memo)
    }
    # refactor the data frame
    x$trv_type <- factor(x$trv_type, levels=mutation_order)

    # sort the data frame so that the duplicated call will remove the
    # proper trv_type
    x <- x[order(x$sample, x$gene, x$trv_type),]

    # collapse the data on sample/gene
    x <- x[!duplicated(x[, c("sample", "gene")]), ]

    return(x)
}
