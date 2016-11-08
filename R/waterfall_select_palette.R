#' waterfall_select_palette
#' 
#' @title Helper function to select a colour palette
#' 
#' @param file_type Which file tyoe is used?
#' @param custom_palette Nullable custom colour palette
#' 
#' @return A vector of colours to be used as palette
waterfall_select_palette <- function(file_type, custom_palette = NULL) {
    if (is.null(custom_palette)) {
        if (toupper(file_type) == toupper('MGI')) {
            palette <- c("nonsense"='#4f00A8', "frame_shift_del"='#A80100',
                         "frame_shift_ins"='#CF5A59', "splice_site_del"='#A80079',
                         "splice_site_ins"='#BC2D94', "splice_site"='#CF59AE',
                         "nonstop"='#000000', "in_frame_del"='#006666',
                         "in_frame_ins"='#00A8A8', "missense"='#009933',
                         "splice_region_del"='#ace7b9', "splice_region_ins"='#cdf0d5',
                         "splice_region"='#59CF74', "5_prime_flanking_region"='#002AA8',
                         "3_prime_flanking_region"='#5977CF',
                         "3_prime_untranslated_region"='#F37812',
                         "5_prime_untranslated_region"='#F2B079', "rna"='#888811',
                         "intronic"='#FDF31C', "silent"='#8C8C8C')
        } else if (toupper(file_type) == toupper('MAF')) {
            palette <- c("Nonsense_Mutation"="grey", "Frame_Shift_Ins"='#A80100',
                         "Frame_Shift_Del"='#CF5A59', "In_Frame_Ins"='#A80079',
                         "In_Frame_Del"='#CF59AE', "Nonstop_Mutation"='#000000',
                         "Translation_Start_Site"='#9159CF', "Splice_Site"='#4f00A8',
                         "Missense_Mutation"='#59CF74', "5\'Flank"='#00A8A8',
                         "3\'Flank"='#79F2F2', "5\'UTR"='#006666',
                         "3\'UTR"='#002AA8', "RNA"='#5977CF', "Intron"='#F37812',
                         "IGR"='#F2B079', "Silent"='#888811',
                         "Targeted_Region"='#FDF31C')
        } else if (toupper(file_type) == toupper('Custom')) {
            memo <- paste0("Defining a palette in mainPalette is recommended ",
                           "when file_type is set to \"Custom\", defaulting to ",
                           "a predefined palette with 20 levels")
            warning(memo)
            palette <- c('#4f00A8', '#A80100', '#CF5A59', '#A80079', '#BC2D94',
                         '#CF59AE', '#000000', '#006666', '#00A8A8', '#009933',
                         '#ace7b9', '#cdf0d5', '#59CF74', '#002AA8', '#5977CF',
                         '#F37812', '#F2B079', '#888811', '#FDF31C', '#8C8C8C')
        }
    }
    else {
        if(toupper(file_type) == "CUSTOM") {
            palette <- custom_palette
        } else {
            if(toupper(file_type) == "MGI") {
                breaks <- c("nonsense", "frame_shift_del",
                            "frame_shift_ins", "splice_site_del",
                            "splice_site_ins", "splice_site",
                            "nonstop", "in_frame_del",
                            "in_frame_ins", "missense",
                            "splice_region_del", "splice_region_ins",
                            "splice_region", "5_prime_flanking_region",
                            "3_prime_flanking_region",
                            "3_prime_untranslated_region",
                            "5_prime_untranslated_region", "rna",
                            "intronic", "silent")
            } else if (toupper(file_type) == "MAF") {
                breaks <- c("Nonsense_Mutation", "Frame_Shift_Ins",
                    "Frame_Shift_Del", "In_Frame_Ins",
                    "In_Frame_Del", "Nonstop_Mutation",
                    "Translation_Start_Site", 
                    "Splice_Site", "Missense_Mutation", 
                    "5\'Flank", "3\'Flank", "5\'UTR",
                    "3\'UTR", "RNA", "Intron",
                    "IGR", "Silent", "Targeted_Region")
            }
            ## Ensure custom palette is correct by
            ## Labelling correctly and making correct length
            palette <- setNames(
                grDevices::colorRampPalette(custom_palette)(length(breaks)), 
                breaks) 
        }
    }
    palette
}