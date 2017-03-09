#' align plots
#'
#' align mutation landscape, mutation burden on sample, and mutation burden on
#' gene plots
#' @name waterfall_align
#' @param heatmap ggplot object displaying a mutation landscape
#' @param genes ggplot object displaying mutation burden on gene
#' @param burden ggplot object displaying mutation burden on sample
#' @param clinical ggplot object displaying clinical information "optional"
#' @param proportions ggplot object displaying proportion of mutation types "optional"
#' @param section_heights Heights of each section (should sum to one)
#' @return a grob object
#' @importFrom gridExtra arrangeGrob

waterfall_align <- function(genes, heatmap, burden, clinical, proportions, 
    section_heights) {
    if (missing(section_heights)) {    
        if (missing(clinical)) {
            if (is.null(proportions)) {
                section_heights <- c(1, 4)
            } else {
                section_heights <- c(1, 4, 0.5)
            }
        } else {
            if (is.null(proportions)) {
                section_heights <- c(1, 4, 0.5)
            } else {
                section_heights <- c(1, 4, 1.2, 0.5)
            }
        }
    }

    # define the ggplot's as grobs and create a blank plot
    genes_grob <- suppressWarnings(ggplot2::ggplotGrob(genes))

    heatmap_grob <- ggplot2::ggplotGrob(heatmap)
    ## Strip out legends and plot separately
    ind_legend <- grep("guide", heatmap_grob$layout$name)
    heatmap_legend <- heatmap_grob[["grobs"]][[ind_legend]]
    heatmap_width <- sum(heatmap_legend$width)
    heatmap_grob[["grobs"]][[ind_legend]] <- zeroGrob()

    burden_grob <- ggplot2::ggplotGrob(burden)
    ## Strip out legends and plot separately
    ind_legend <- grep("guide", burden_grob$layout$name)
    burden_legend <- burden_grob[["grobs"]][[ind_legend]]
    burden_width <- sum(burden_legend$width)
    burden_grob[["grobs"]][[ind_legend]] <- zeroGrob()

    ## Strip out y axis titles and plot separately
    ind_axis <- grep("ylab-l", burden_grob$layout$name)
    burden_axis_title <- burden_grob[["grobs"]][[ind_axis]]
    # burden_axis_title$children[[1]]$vjust <- 1
    print(burden_axis_title$children[[1]]$vjust)
    burden_axis_title$children[[1]]$hjust <- 0.5
    burden_grob[["grobs"]][[ind_axis]] <- zeroGrob()

    if (!is.null(proportions)) {
        prop_grob <- ggplot2::ggplotGrob(proportions)
        ind_legend <- grep("guide", prop_grob$layout$name)
        prop_legend <- prop_grob[["grobs"]][[ind_legend]]
        prop_width <- sum(prop_legend$width)
        prop_grob[["grobs"]][[ind_legend]] <- zeroGrob()
        
        ## Strip out y axis titles and plot separately
        ind_axis <- grep("ylab-l", prop_grob$layout$name)
        prop_axis_title <- prop_grob[["grobs"]][[ind_axis]]
        print(prop_axis_title$children[[1]]$vjust)
        prop_axis_title$children[[1]]$vjust <- 0
        prop_axis_title$children[[1]]$hjust <- -0.2
        prop_grob[["grobs"]][[ind_axis]] <- zeroGrob()

    } else prop_width <- NULL
    blankPanel <- grid::grid.rect(gp=grid::gpar(col="white"))
    if(!missing(clinical)) {
        clin_grob <- ggplot2::ggplotGrob(clinical)
        ind_legend <- grep("guide", clin_grob$layout$name)
        clin_legend <- clin_grob[["grobs"]][[ind_legend]]
        clin_width <- sum(clin_legend$width)
        clin_grob[["grobs"]][[ind_legend]] <- zeroGrob()

        ## Strip out y axis titles and plot separately
        ind_axis <- grep("ylab-l", clin_grob$layout$name)
        clin_axis_title <- clin_grob[["grobs"]][[ind_axis]]
        clin_grob[["grobs"]][[ind_axis]] <- zeroGrob()

    } else clin_width <- NULL

    legend_width <- unit.pmax(heatmap_width, burden_width, clin_width, prop_width)
    width_left <- unit(1, "npc") - max(legend_width)
    widths <- grid::unit.c(width_left * 0.15, width_left * 0.85, legend_width)

    # Adjust the grob widths so heatmap and burden plots line up
    if(!missing(clinical))
    {
        if (!is.null(proportions)) {
            maxwidth <- grid::unit.pmin(heatmap_grob$widths,
                                       burden_grob$widths,
                                       clin_grob$widths,
                                       prop_grob$widths)
            prop_grob$widths <- as.list(maxwidth)
        } else {
            maxwidth <- grid::unit.pmin(heatmap_grob$widths,
                                       burden_grob$widths,
                                       clin_grob$widths
                                       )
        }
        burden_grob$widths <- as.list(maxwidth)
        heatmap_grob$widths <- as.list(maxwidth)
        clin_grob$widths <- as.list(maxwidth)
    } else {
        if (!is.null(proportions)) {
            maxwidth <- grid::unit.pmin(heatmap_grob$widths,
                                       burden_grob$widths,
                                       prop_grob$widths)
            prop_grob$widths <- as.list(maxwidth)
        } else {   
            maxwidth <- grid::unit.pmin(heatmap_grob$widths,
                                       burden_grob$widths
                                       )
        }
        burden_grob$widths <- as.list(maxwidth)
        heatmap_grob$widths <- as.list(maxwidth)
    }

    # Adjust the grob heights so heatmap, and genes plots line up
    maxheight <- grid::unit.pmax(genes_grob$heights, heatmap_grob$heights)
    genes_grob$heights <- as.list(maxheight)
    heatmap_grob$heights <- as.list(maxheight)

    # plot the grobs with grid.arrange
    if(!missing(clinical)) {
        if (!is.null(proportions)) {
            nrows <- 4
            grobs <- list(
                burden_axis_title, burden_grob, burden_legend, 
                genes_grob, heatmap_grob, heatmap_legend, 
                prop_axis_title, prop_grob, prop_legend, 
                clin_axis_title, clin_grob, clin_legend)
        } else {
            nrows <- 3
            grobs <- list(
                burden_axis_title, burden_grob, burden_legend, 
                genes_grob, heatmap_grob, heatmap_legend,
                clin_axis_title, clin_grob, clin_legend)
        }
        heatmap <- gridExtra::arrangeGrob(grobs = grobs,
                                     ncol=3, nrow=nrows, 
                                     widths=widths,
                                     heights=section_heights)
        
        
        
    } else {
        if (!is.null(proportions)) {
            nrows <- 3
            grobs <- list(
                burden_axis_title, burden_grob, burden_legend,
                genes_grob, heatmap_grob, heatmap_legend, 
                prop_axis_title, prop_grob, prop_legend)
        } else {
            grobs <- list(
                burden_axis_title, burden_grob, burden_legend,
                genes_grob, heatmap_grob, heatmap_legend,
                )
            nrows <- 2
        }
        heatmap <- gridExtra::arrangeGrob(grobs = grobs,
                                     ncol=3, nrow=nrows,
                                     widths=widths, heights=section_heights)
    }

    return(heatmap)
}
