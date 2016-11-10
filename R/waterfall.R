#' Construct a waterfall plot
#'
#' Given a data frame construct a water fall plot showing the mutation burden
#' and mutation type on a gene and sample level.
#' @name waterfall
#' @param x Object of class data frame representing annotated mutations. The
#' data frame supplied must have one of the following sets of column names
#' ("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification") for
#' fileType="MAF", ("sample","gene_name","trv_type") for fileType="MGI" or
#' ("sample", "gene", "variant_class") for fileType="Custom". This columns
#' should represent samples in a cohort, gene with mutation, and the mutation
#' type respectively.
#' @param mainRecurCutoff Numeric value between 0 and 1 specifying a
#' mutation recurrence cutoff. Genes which do not have mutations in the
#' proportion os samples defined are removed.
#' @param mainGrid Boolean specifying if a grid should be overlayed on the main
#' plot. Not recommended if the number of genes or samples to be plotted
#' is large.
#' @param mainXlabel Boolean specifying whether to label the x-axis with sample
#' names. Not recommended if the number of samples to be plotted is large.
#' @param main_geneLabSize Intenger specifying the size of gene names displayed
#' on the y-axis.
#' @param mainLabelCol Character string specifying a column name from the
#' argument supplied to parameter `x` from which to derive cell labels from
#' (see details and vignette).
#' @param mainLabelSize Integer specifying the size of text labels for cells
#' in the main plot. Valid only if argument is supplied to the parameter
#' `mainLabelCol`.
#' @param mainLabelAngle Integer specifying the degree of rotation for
#' text labels. Valid only if argument is supplied to the parameter
#' `mainLabelCol`.
#' @param mainDropMut Boolean specifying whether to drop unused
#' "mutation type" levels from the legend.
#' @param mainPalette Character vector specifying colours for mutation types
#' plotted in the main plot, must specify a colour for each mutation type
#' plotted.
#' @param mainLayer Valid ggplot2 layer to be added to the main plot.
#' @param mutBurden Object of class data frame containing columns "sample",
#' "mut_burden" with sample levels matching those supplied in x.
#' @param plotMutBurden Boolean specify if the mutation burden sub-plot should
#' be displayed.
#' @param coverageSpace Integer specifying the size in bp of the genome
#' covered by sequence data from which mutations could be called
#' (see details and vignette).
#' @param mutBurdenLayer Valid ggplot2 layer to be added to the top sub-plot.
#' @param clinData Object of class data frame with rows representing clinical
#' data. The data frame should be in "long format" and columns must be names as
#' "sample", "variable", and "value" (optional see details and vignette).
#' @param clinLegCol Integer specifying the number of columns in the legend for
#' the clinical data, only valid if argument is supplied to parameter clinData.
#' @param clinVarOrder Character vector specifying the order in which to plot
#' variables in the variable column of the argument given to the parameter
#' clinData. The argument supplied to this parameter should have the same unique
#' length and values as in the variable column of the argument supplied to
#' parameter clinData (see vignette).
#' @param clinVarCol Named character vector specifying the mapping of colours
#' to variables in the variable column of the data frame supplied to clinData
#' (ex. "variable"="colour").
#' @param clinLayer Valid ggplot2 layer to be added to the clinical sub-plot.
#' @param sampRecurLayer Valid ggplot2 layer to be added to the left sub-plot.
#' @param plotGenes Character vector specifying genes to plot. If not null genes
#' not specified within this character vector are removed.
#' @param geneOrder Character vector specifying the order in which to plot
#' genes.
#' @param plotSamples Character vector specifying samples to plot. If not null
#' all other samples not specified within this parameter are removed.
#' @param sampOrder Character vector specifying the order of the samples to
#' plot.
#' @param maxGenes Integer specifying the maximum number of genes to be plotted.
#' Genes kept will be choosen based on the reccurence of mutations in samples.
#' @param rmvSilent Boolean specifying if silent mutations should be removed
#' from the plot.
#' @param fileType Character string specifying the file format of the data
#' frame specified to parameter `x`, one of "MGI", "MAF", "Custom"
#' (see details and vignette).
#' @param variant_class_order Character vector specifying the hierarchical order
#' of mutation types to plot, required if file_type == "Custom"
#' (see details and vignette).
#' @param out Character vector specifying the the object to output, one of
#' "data", "grob", or "plot", defaults to "plot" (see returns).
#' @param plot_proportions Plot mutational profile layer?
#' @param proportions_layer ggplot2 layer(s) to be added to the mutational 
#'      profile plot
#' @param proportions_type Which type of proportions plot to use? Can be 
#'  "trv_type" or "TvTi" currently
#' @param section_heights Heights of each section. Must be the same length as 
#'  the number of vertical section
#' @details waterfall is a function designed to visualize the mutations seen in
#' a cohort. The function takes a data frame with appropriate column names (see
#' fileType parameter) and plots the mutations within. In cases where multiple
#' mutations occur in the same cell the most deleterious mutation is given
#' priority (see vignette for default priority). If the fileType parameter is
#' set to "Custom" the user most supply this priority via the
#' `variant_class_order` parameter with the highest priorities occuring first.
#' Additionally this parameter will override the default orders of MGI and MAF
#' file types.
#'
#' Various data subsets are allowed via the waterfall function (see above), all
#' of these subsets will occur independently of the mutation burden calculation.
#' To clarify the removal of genes and mutations will only occur after the
#' mutation burden is calculated. The mutation burden calculation is only meant
#' to provide a rough estimate and assumes that the coverage breadth within the
#' cohort is aproximately equal. For more accurate calculations it is
#' recommended to supply this information via the mutBurden parameter which.
#' Note that the mutation burden calculation relies on the `coverageSpace`
#' parameter (see vignette).
#'
#' It is possible to display additional information within the plot via cell
#' labels. The `mainLabelCol` parameter will look for an additional column in
#' the data frame and plot text within cells based on those values
#' (see vignette).
#' @examples
#' # Plot the data
#' waterfall(brcaMAF, plotGenes=c("PIK3CA", "TP53", "USH2A", "MLL3", "BRCA1"))
#' @return One of the following, a list of dataframes containing data to be
#' plotted, a grob object, or a plot.
#' @importFrom utils tail
#' @export

waterfall <- function(x, mainRecurCutoff=0, mainGrid=TRUE, mainXlabel=FALSE, 
                      main_geneLabSize=8, mainLabelCol=NULL, mainLabelSize=4,
                      mainLabelAngle=0, mainDropMut=FALSE, mainPalette=NULL, 
                      mainLayer=NULL, mutBurden=NULL, plotMutBurden=TRUE, 
                      coverageSpace=44100000, mutBurdenLayer=NULL,
                      clinData=NULL, clinLegCol=1, clinVarOrder=NULL,
                      clinVarCol=NULL, clinLayer=NULL, sampRecurLayer=NULL, 
                      plotGenes=NULL, geneOrder=NULL, plotSamples=NULL,
                      sampOrder=NULL, maxGenes=NULL, rmvSilent=FALSE,
                      fileType='MAF', variant_class_order=NULL, out="plot",
                      plot_proportions = FALSE,
                      proportions_layer = NULL, proportions_type = "TRV_TYPE",
                      section_heights)
{
  
    # Perform data quality checks and conversions
    inputDat <- waterfall_qual(x, clinData, mutBurden, file_type=fileType,
                               label_col=mainLabelCol)
    data_frame <- inputDat[[1]]
    clinData <- inputDat[[2]]
    mutBurden <- inputDat[[3]]

    # Set flag if it is desirable to plot cell text
    if(!is.null(mainLabelCol))
    {
        main.plot_label_flag <- TRUE
    } else {
        main.plot_label_flag <- FALSE
    }

    # If it is requested subset the input data on a sample list
    if(!is.null(plotSamples))
    {
        data_frame <- waterfall_sampAlt(data_frame, plotSamples)
    }

    # add in a count of mutations at the sample level before anything is
    # stripped out and save for mutation recurrence plot
    data_frame2 <- waterfall_calcMutFreq(data_frame[,c('sample', 'gene',
                                                       'trv_type')])

    # Subset the data to remove silent mutations if specified
    if(rmvSilent==TRUE)
    {
        data_frame <- waterfall_rmvSilent(data_frame)
    }

    # Subset the data based on a vector of genes if supplied
    if(!is.null(plotGenes))
    {
        data_frame <- waterfall_geneAlt(data_frame, plotGenes)
    }

    # Remove trv_type that are not the most deleterious for a given gene/sample
    data_frame <- waterfall_hierarchyTRV(data_frame, fileType,
                                         variant_class_order)

    # Subset the data based on the recurrence of mutations at the gene level
    data_frame <- waterfall_geneRecurCutoff(data_frame, mainRecurCutoff)

    # Use the max genes parameter to limit the number of genes plotted
    # and then reorder genes based on the frequency of mutations
    gene_sorted <- waterfall_geneSort(data_frame, geneOrder)
    data_frame$gene <- factor(data_frame$gene, levels=gene_sorted)
    if(!is.null(maxGenes))
    {
        max_gene_list <- utils::tail(gene_sorted, maxGenes)
        data_frame <- data_frame[data_frame$gene %in% max_gene_list,]
    }

    # reorder the samples based on hiearchial sort on ordered gene list
    sample_order <- waterfall_sampSort(data_frame, sampOrder)
    data_frame$sample <- factor(data_frame$sample, levels=sample_order)

    if (plot_proportions) {
        prop_x_label <- is.null(clinData)
        if (toupper(proportions_type) == "TRV_TYPE") {
            prop_dat <- inputDat[[1]]
            prop_dat[["sample"]] <- factor(prop_dat[["sample"]], 
              levels = sample_order)
            p5 <- waterfall_build_proportions(
                data_frame = prop_dat, 
                plot_palette = mainPalette,
                file_type = fileType,
                layers = proportions_layer,
                x_label = prop_x_label
            )
        } else if (toupper(proportions_type) == "TVTI") {
            if (is.null(proportions_layer)) {
                proportions_layer <- list(
                    theme(
                        axis.ticks.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.text.y = element_blank(),
                        panel.background = element_blank(),
                        panel.border = element_blank(),
                        panel.grid.minor = element_blank(),
                        plot.background = element_blank()
                    ),
                    scale_y_continuous(
                        name = "Proportion", 
                        labels = scales::percent_format()
                    ),
                    guides(fill = guide_legend(ncol = 2))
                )
            }
            p5 <- TvTi(
                x,
                fileType = fileType,
                sort = "custom",
                sample_order_input = sample_order,
                layers = proportions_layer,
                return_plot = TRUE
            )
            if (prop_x_label) p5 <- p5 + xlab(paste0('Sample (n=', nlevels(data_frame$sample), ')'))
        }
    } else p5 <- NULL

    # Reorder the sample levels in data_frame2 to match the main plot's levels,
    # and then plot the top margin plot
    if(isTRUE(plotMutBurden))
    {
        if(!is.null(mutBurden))
        {
            if(!setequal(sample_order, mutBurden$sample))
            {
                stop("levels in the sample column of mutBurden does not match
                 either: the samples given in x, or plotSamples")
            }

            mutBurden$sample <- factor(mutBurden$sample, levels=sample_order)
            p3 <- waterfall_buildMutBurden_B(mutBurden, layers=mutBurdenLayer)
        } else {
            data_frame2$sample <- factor(data_frame2$sample,
                                         levels=sample_order)
            p3 <- waterfall_buildMutBurden_A(data_frame2, coverageSpace,
                                             layers=mutBurdenLayer)
        }
    } else {
        # create a blank ggplot object
        df <- data.frame()
        p3 <- ggplot2::ggplot(df) + ggplot2::geom_point() +
            ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
            ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                           axis.text.y=ggplot2::element_blank(),
                           axis.ticks.x=ggplot2::element_blank(),
                           axis.ticks.y=ggplot2::element_blank(),
                           panel.background=ggplot2::element_blank(),
                           panel.grid=ggplot2::element_blank())
    }

    # Plot the Left Bar Chart
    p2 <- waterfall_buildGenePrevelance(data_frame, layers=sampRecurLayer)

    # if there are any NA values in the data frame for a gene, give these NA
    # values a gene name so they are plotted properly
    data_frame <- waterfall_NA2gene(data_frame)

    # Plot the Heatmap
    plot_x_title <- !(!is.null(clinData) || plot_proportions)

    p1 <- waterfall_buildMain(data_frame, grid=mainGrid,
                              label_x=mainXlabel,
                              gene_label_size=main_geneLabSize,
                              file_type=fileType,
                              drop_mutation=mainDropMut,
                              plot_x_title=plot_x_title,
                              plot_label=main.plot_label_flag,
                              plot_label_size=mainLabelSize,
                              plot_palette=mainPalette, layers=mainLayer,
                              plot_label_angle=mainLabelAngle)
    

    # Plot any clinical data if it is specified
    if(!is.null(clinData))
    {
        # match the levels of sample in y to conform to the main plot
        clinData$sample <- factor(clinData$sample, levels=sample_order)
        if(any(is.na(clinData$sample)))
        {
            clinData <- clinData[-which(is.na(clinData$sample)),]
        }

        # plot the clinical data
        p4 <- multi_buildClin(clinData, clin.legend.col=clinLegCol,
                              clin.var.colour=clinVarCol,
                              clin.var.order=clinVarOrder,
                              clin.layers=clinLayer)

        # Align all plots and return as 1 plot
        pA <- waterfall_align(p2 = p2, p1 = p1, p3 = p3, p4 = p4, p5 = p5, 
          section_heights = section_heights)
        return(grid::grid.draw(pA))
    } else 
    {
        pA <- waterfall_align(p2 = p2, p1 = p1, p3 = p3, p5 = p5,
          section_heights = section_heights)
    }
    dataOut <- list("main"=data_frame,
                    "mutation_count"=data_frame2,
                    "clinical"=clinData)
    output <- multi_selectOut(data=dataOut, plot=pA, draw=TRUE, out=out)
    return(output)
}
