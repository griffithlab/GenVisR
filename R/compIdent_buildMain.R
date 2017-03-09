#' Compare sample identities
#'
#' Produce an identity SNP plot displaying VAFs of 24 SNP locations and coverage
#' information to compare multiple sample identities
#' @name compIdent_buildMain
#' @param x Data frame of vaf for each sample
#' @param mainLayer Valid ggplot2 layer for altering the main plot.
#' @param covLayer Valid ggplot2 layer for altering the coverage plot.
#' @return ggplot2 grob object
#' @import ggplot2
#' @importFrom grid unit.pmax
#' @importFrom gridExtra rbind.gtable
#' @importFrom gridExtra arrangeGrob
#' @importFrom gtable gtable_add_rows
#' @importFrom gtable gtable_add_grob
#' @importFrom scales log10_trans

compIdent_buildMain <- function(x, mainLayer=NULL, covLayer=NULL)
{
    ############### Plot VAF for main plot ###################
    varBases <- x[x$sample == x$sample[1],]$var
    refBases <- x[x$sample == x$sample[1],]$getref
    
    main <- ggplot(x, aes_string(x='name', y='vaf', 
                                 colour='sample', fill='sample'))
    plotVAF <- geom_point(size=5, position=position_dodge(width=0.5))
    y_axis <- scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1))
    y_label <- ylab("VAF")
    varBases <- x[x$sample == x$sample[1],]$var
    x_axis <- scale_x_discrete(labels=paste0(refBases, "/", varBases),
                               position="top", name="Reference/Variant")
    plot_theme_A <- theme_bw()
    plot_theme_B <- theme(legend.position="none",
                          panel.grid.minor.y=element_blank(),
                          panel.grid.minor.x=element_line(colour="grey"),
                          plot.margin=unit(c(1, 1, 0, 1), "cm"),
                          axis.ticks=element_blank())
    
    # set layer and grab the size of the axis.text.x theme param if it exists
    if(!is.null(mainLayer)){
        layerA <- mainLayer
        x_axis_text_size <- layerA$axis.text.x$size
    } else {
        layerA <- geom_blank()
    }

    p1 <- main + plotVAF + y_axis + x_axis + y_label + plot_theme_A +
        plot_theme_B + layerA
    
    ############# produce a coverage plot ##################
    main2 <- ggplot(x, aes_string(x='name', y='total_reads',
                                  fill='sample', colour='sample'))
    plotCov <- geom_bar(stat="identity", position="dodge")
    x_label <- xlab("SNP")
    y_label <- ylab("Coverage")
    y_axis <- scale_y_continuous(trans=scales::log10_trans(),
                                breaks=c(10, 50, 100, 500, 1000))
    theme_A <- theme(legend.position="bottom", axis.text.y=element_text(hjust=0.5),
                     axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5),
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x=element_line(colour="grey"),
                     plot.margin=unit(c(0,1,1,1), "cm"))
    
    # set layer
    if(!is.null(covLayer)){
        layerB <- covLayer
    } else {
        layerB <- geom_blank()
    }
    
    p2 <- main2 + plotCov + x_label + y_label + y_axis + theme_bw() +
        theme_A + layerB

    ############### Ensure plot widths are identical and plot ################
    p2 <- ggplotGrob(p2)
    p1 <- ggplotGrob(p1)
    maxwidth = grid::unit.pmax(p1$widths, p2$widths)
    p1$widths <- as.list(maxwidth)
    p2$widths <- as.list(maxwidth)

    # Plot all
    final <- gridExtra::arrangeGrob(p1, p2, heights=c(60, 40))

    return(final)
}
