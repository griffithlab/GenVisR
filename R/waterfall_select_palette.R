#' select_palette
#' 
#' @param file_type
#' 
#' @return A vector of colours to be used as palette
select_palette <- function(file_type) {
    if(toupper(file_type) == toupper('MGI')) {
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
    } else if(toupper(file_type) == toupper('MAF')) {
        palette <- c("Nonsense_Mutation"="grey", "Frame_Shift_Ins"='#A80100',
                     "Frame_Shift_Del"='#CF5A59', "In_Frame_Ins"='#A80079',
                     "In_Frame_Del"='#CF59AE', "Nonstop_Mutation"='#000000',
                     "Translation_Start_Site"='#9159CF', "Splice_Site"='#4f00A8',
                     "Missense_Mutation"='#59CF74', "5\'Flank"='#00A8A8',
                     "3\'Flank"='#79F2F2', "5\'UTR"='#006666',
                     "3\'UTR"='#002AA8', "RNA"='#5977CF', "Intron"='#F37812',
                     "IGR"='#F2B079', "Silent"='#888811',
                     "Targeted_Region"='#FDF31C')
    } else if(toupper(file_type) == toupper('Custom')) {
        memo <- paste0("Defining a palette in mainPallete is recommended ",
                       "when file_type is set to \"Custom\", defaulting to ",
                       "a predefined palette with 20 levels")
        warning(memo)
        palette <- c('#4f00A8', '#A80100', '#CF5A59', '#A80079', '#BC2D94',
                     '#CF59AE', '#000000', '#006666', '#00A8A8', '#009933',
                     '#ace7b9', '#cdf0d5', '#59CF74', '#002AA8', '#5977CF',
                     '#F37812', '#F2B079', '#888811', '#FDF31C', '#8C8C8C')
    }
    palette
}