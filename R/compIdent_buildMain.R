#' Compare sample identities
#'
#' Produce an identity SNP plot displaying VAFs of 24 SNP locations and coverage
#' information to compare multiple sample identities
#' @name compIdent_buildMain
#' @param count_tables list of readcount tables
#' @param sample_names list of names associated with each sample in counts list
#' @return ggplot2 object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom gtable gtable_add_rows
#' @importFrom gridExtra arrangeGrob

compIdent_buildMain <- function(count_tables, sample_names)
{
    # Function that calculates frequencies from readcounts.
    freq <- function(v){
        v[,5:8] <- v[,5:8]/v[,4]
        return(v)
    }
    freqs <- lapply(count_tables, freq)

    # Add a column of the variant allele after the column of the
    # reference allele
    var <- c('A','A','C','T','G','C','C','A','C','C','G','A','G','T','G','A','A',
             'A','A','C','C','G','G','C')
    add_varcol <- function(w)
    {
        cbind(w[,1:3],var,w[,4:8])
    }
    freqs_withvarcol <- lapply(freqs, add_varcol)

    # Extract VAF based upon variant allele
    calcVAF <- function(x)
    {
        apply(x,1,function(y){
              return(as.numeric(y[as.character(y["var"])]))
              })
    }
    calcVAFs <- lapply(freqs_withvarcol, calcVAF)

    # Make one table: rsid, sampleVAFs
    id <- 1:24
    vaf <- cbind(id, as.data.frame(calcVAFs))
    vaf.long <- reshape2::melt(vaf, id.vars='id')

    # factor variable column
    vaf.long$variable <- factor(vaf.long$variable)
    
    # Plot VAF
    # The point represents the VAF (y-axis) for each dbSNP (x-axis). The colour of
    # the points represent samples. Title "Variant Base" refers to the x-axis
    # labels when they are moved above the graph.
    plotVAF <- ggplot(vaf.long, aes_string(x='id', y='value')) +
    geom_point(aes_string(colour='variable'), size=5,
               position=position_dodge(width=0.5)) +
    ggtitle("Variant Base")
    
    # y-axis: label as "VAF". Set major breaks at 0, 0.5, and 1.0.
    plotVAF <- plotVAF + ylab("VAF") + scale_y_continuous(limits=c(0,1),
                                                          breaks=c(0,0.50,1.00))
    
    # x-axis: Show the variant allele at each locus as x-axis tick labels.
    plotVAF <- plotVAF + scale_x_continuous(limits=c(0.5,24.5), breaks=1:24,
                                            labels=as.list(var))
    
    # theme: Change panel to white. Remove the legend
    # Remove the minor y-axis lines. Add grey lines between each dbSNP to separate
    # Remove x-axis ticks
    plotVAF <- plotVAF + theme_bw() +
    theme(legend.position="none", panel.grid.minor.y=element_blank(),
          panel.grid.minor.x=element_line(colour="grey"),
          plot.margin=unit(c(1,1,0,1), "cm"),
          axis.title.x=element_blank(),
          axis.ticks=element_blank())

    ## Move x-axis labels above graph
    plotVAF <- ggplotGrob(plotVAF)
    panel <- c(subset(plotVAF$layout, `name`=="panel", se=t:r))    
    rn <- which(plotVAF$layout$name == "axis-b")
    axis.grob <- plotVAF$grobs[[rn]]
    axisb <- axis.grob$children[[2]]
    axisb
    axisb$heights <- rev(axisb$heights)
    axisb$grobs <- rev(axisb$grobs)
    plotVAF <- gtable::gtable_add_rows(plotVAF,
                                       plotVAF$heights[plotVAF$layout[rn, ]$l],
                                       panel$t-1)
    plotVAF <- gtable::gtable_add_grob(plotVAF, axisb, l=panel$l, t=panel$t,
                                       r=panel$r)
    plotVAF <- plotVAF[-(panel$b+2),]

    ## Add Reference Base underneath plotVAF

    # Data frame (1 row by 24 columns) of the reference bases at each dbSNP
    a <- rep(1,24)
    
    # Extract reference alleles from one of the counts table (column 3)
    ref <- count_tables[[1]][,3]
    .df <- data.frame(vaf[1], ref, a)

    # Create a blank plot for place-holding
    plot.df <- ggplot(.df, aes(x=id,y=a, label=ref)) + geom_blank() +
    xlab("Reference Base") +
    scale_x_continuous(limits=c(0.5,24.5), breaks=1:24) +
    scale_y_continuous(limits=c(0.99,1.01)) +
    theme_bw() +
    theme(plot.margin=unit(c(0,1,0,1), "cm"),
          axis.title.x=element_text(size=14, vjust=1),
          axis.title.y=element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank())

    plot.df <- plot.df + geom_text(size=3.5)

    # Create a data frame containing only total readcounts
    totalcounts <- function(y){
        apply(y,1,function(z){
              return(as.vector(z["total_reads"]))
              })
    }
    readcounts <- lapply(count_tables, totalcounts)
    readcounts.df <- cbind(id, as.data.frame(readcounts))
    
    ## Use colnames to specify sample names, which will appear in legend
    colnames(readcounts.df) <- c('id', sample_names)
    readcounts.long <- reshape2::melt(readcounts.df, id.vars='id',
                                      measure.vars=sample_names)

    # Plot read counts
    dbSNP_rsID <- c('rs2229546','rs1410592','rs497692','rs10203363','rs2819561',
                    'rs4688963','rs309557','rs2942','rs17548783','rs4735258',
                    'rs1381532','rs10883099','rs4617548','rs7300444',
                    'rs9532292','rs2297995','rs4577050','rs2070203','rs1037256',
                    'rs9962023','rs2228611','rs10373','rs4148973','rs4675')

    readcounts.long$value <- as.numeric(readcounts.long$value)
    plotRC <- ggplot(readcounts.long,
                     aes_string(x='id', y='value', fill='variable')) +
    geom_bar(stat="identity",position="dodge")

    # x-axis: Label as "dbSNP rsID". Label ticks with rsIDs.
    plotRC <- plotRC + xlab("dbSNP rsID") + scale_x_continuous(limits=c(0.5,24.5),
                                                               breaks=1:24,
                                                               labels=dbSNP_rsID)

    # y-axis: Label as "Total Readcount". Make a logarithmic scale. Set breaks.
    plotRC <- plotRC + ylab("Total Readcount") +
    scale_y_continuous(trans=log10_trans(),
                       breaks=c(10,50,100,500,1000))

    # theme: Make panel white, Move legend underneath graph, Rotate x-axis labels
    # 90deg
    plotRC <- plotRC + theme_bw() +
    theme(legend.position="bottom", axis.text.y=element_text(hjust=0.5),
          axis.text.x=element_text(angle=90, vjust=0.5, hjust=0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x=element_line(colour="grey"),
          plot.margin=unit(c(0,1,1,1), "cm"))

    ## Graph plotVAF over plot.df over plotRC
    # Change widths of plot.df and plotRC to width of plotVAF
    plot.df <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot.df))
    plot.df$widths <- plotVAF$widths
    plotRC <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plotRC))
    plotRC$widths <- plotVAF$widths

    # Plot all
    plotVAF_refbases <- gridExtra::arrangeGrob(plotVAF, plot.df, heights=c(88,12))
    plotVAFoverRC <- gridExtra::arrangeGrob(plotVAF_refbases, plotRC,
                                            heights=c(60,40))

    return(grid::grid.draw(plotVAFoverRC))
}
