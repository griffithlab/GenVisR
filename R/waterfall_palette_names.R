#'  waterfall_palette_names
#' 
#' @title waterfall_palette_names
#' 
#' @description Make labels and breaks for palettes
#' 
#' @param palette Named colour vector as input 
#' @param file_type Which file type is involved?
#' @param data_frame Only used if file_type is "custom"
#' 
#' @return a named list of "breaks" and "labels"
waterfall_palette_names <- function(palette, file_type, data_frame) {
    # Create breaks specific and labels for specified file type
    ## Create labels and breaks from 
    ## names of palette to avoid ordering issues       
    if(toupper(file_type) == toupper('MGI'))
    {
        breaks <- names(palette)

        labels <- gsub("(\\d)_prime", "\\1'", breaks)
        labels <- gsub("untranslated_region", "UTR", labels)
        labels <- gsub("flanking_region", "Flank", labels)
        labels <- gsub("_", " ", labels)
        labels <- gsub("del$", "deletion", labels)
        labels <- gsub("ins$", "insertion", labels)
        labels <- gsub("prime ", "", labels)
        labels <- gsub("(rna)", "\\U\\1", labels, perl = TRUE)
        labels <- gsub("\\b([a-z])", "\\U\\1", labels, perl = TRUE)
    } else if(toupper(file_type) == toupper('MAF')) {
        breaks <- names(palette)
        labels <- gsub("_", " ", breaks)
        labels <- gsub("'", "' ", labels)
    } else if(toupper(file_type) == toupper('Custom')) {
        breaks <- levels(data_frame[["trv_type"]])
        labels <- breaks
    }
    list("breaks" = breaks,
        "labels" = labels)
}

