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
#' @importFrom grid nullGrob
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
    
    message("This function has been deprecated in order to implement an object oriented programming style! Please use Waterfall() with a capital W instead!")
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
            proportions_plot <- waterfall_build_proportions(
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
            proportions_plot <- TvTi(
                x,
                fileType = fileType,
                sort = "custom",
                sample_order_input = sample_order,
                layers = proportions_layer,
                return_plot = TRUE
            )
            if (prop_x_label) proportions_plot <- proportions_plot + xlab(paste0('Sample (n=', nlevels(data_frame$sample), ')'))
        }
    } else proportions_plot <- NULL
    
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
            burden_plot <- waterfall_buildMutBurden_B(mutBurden, layers=mutBurdenLayer)
            } else {
                data_frame2$sample <- factor(data_frame2$sample,
                                             levels=sample_order)
                burden_plot <- waterfall_buildMutBurden_A(data_frame2, coverageSpace,
                                                          layers=mutBurdenLayer)
            }
    } else {
        # create a blank ggplot object
        burden_plot <- grid::nullGrob()
    }
    
    # Plot the Left Bar Chart
    gene_plot <- waterfall_buildGenePrevelance(data_frame, layers=sampRecurLayer,
                                               gene_label_size=main_geneLabSize)
    
    # if there are any NA values in the data frame for a gene, give these NA
    # values a gene name so they are plotted properly
    data_frame <- waterfall_NA2gene(data_frame)
    
    # Plot the Heatmap
    plot_x_title <- !(!is.null(clinData) || plot_proportions)
    
    heatmap <- waterfall_buildMain(data_frame, grid=mainGrid,
                                   label_x=mainXlabel,
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
        clinical_plot <- multi_buildClin(clinData, clin.legend.col=clinLegCol,
                                         clin.var.colour=clinVarCol,
                                         clin.var.order=clinVarOrder,
                                         clin.layers=clinLayer)
        
        # Align all plots and return as 1 plot
        pA <- waterfall_align(genes = gene_plot, heatmap = heatmap, burden = burden_plot, 
                              clinical = clinical_plot, proportion = proportions_plot, 
                              section_heights = section_heights)
        return(grid::grid.draw(pA))
    } else 
    {
        pA <- waterfall_align(genes = gene_plot, heatmap = heatmap, burden = burden_plot, 
                              proportion = proportions_plot, section_heights = section_heights)
    }
    dataOut <- list("main"=data_frame,
                    "mutation_count"=data_frame2,
                    "clinical"=clinData)
    output <- multi_selectOut(data=dataOut, plot=pA, draw=TRUE, out=out)
    return(output)
}
#' Convert Custom File
#'
#' Convert columns of a Custom annotation file into a format
#' recognizable by internal functions
#' @name waterfall_Custom2anno
#' @param x a data frame with columns having values for sample, gene, mutation
#' type
#' @param label_col Character string specifying the column name of a
#' label column (optional)
#' @noRd
#' @return a data frame coerced from custom to annotation format

waterfall_Custom2anno <- function(x, label_col)
{
    # message statement
    memo <- paste0("Detected \"Custom\" file_type flag, ",
                   "looking for correct column names...")
    message(memo)
    
    # define expected columns
    expec_col <- c("sample", "gene", "variant_class")
    if(!is.null(label_col))
    {
        expec_col <- c(expec_col, label_col)
    }
    
    # check expected columns are present
    if(!all(expec_col %in% colnames(x)))
    {
        memo <- paste0("Did not detect correct column names, column names
                       should be: ", toString(expec_col))
        stop(memo)
    }
    
    x <- x[,c('sample', 'gene', 'variant_class', label_col)]
    
    if(!is.null(label_col))
    {
        colnames(x) <- c('sample', 'gene', 'trv_type', 'label')
    } else {
        colnames(x) <- c('sample', 'gene', 'trv_type')
    }
    
    # if no silent mutations are present warn the user
    if(all(!toupper(x$trv_type) %in% toupper("silent")))
    {
        warning("Did not detect silent mutations in input, is this expected?")
    }
    return(x)
}
#' Convert MAF File
#'
#' Convert columns of a mutation annotation file "MAF" into a format
#' recognizable by internal functions
#' @name waterfall_MAF2anno
#' @param x a data frame in MAF format
#' @param label_col Character string specifying the column name of a
#' label column
#' @noRd
#' @return a data frame coerced from MAF to TGI format

waterfall_MAF2anno <- function(x, label_col)
{
    # Check that correct column names are present and convert to internal format
    expec_col <- c('Tumor_Sample_Barcode', 'Hugo_Symbol',
                   'Variant_Classification')
    
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
    
    x <- x[,c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification',
              label_col)]
    
    if(!is.null(label_col))
    {
        colnames(x) <- c('sample', 'gene', 'trv_type', 'label')
    } else {
        colnames(x) <- c('sample', 'gene', 'trv_type')
    }
    return(x)
}
#' Convert MGI File
#'
#' Convert columns of a mutation annotation file "MGI" into a format
#' recognizable by internal functions
#' @name waterfall_MGI2anno
#' @param x a data frame in MGI internal format
#' @param label_col Character string specifying the column name of a label
#' column
#' @noRd
#' @return a data frame coerced from MGI to internal annotation format

waterfall_MGI2anno <- function(x, label_col)
{
    # Check that correct column names are present and convert to internal format
    expec_col <- c('sample', 'gene_name', 'trv_type')
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
    
    x <- x[,c('sample', 'gene_name', 'trv_type', label_col)]
    if(!is.null(label_col))
    {
        colnames(x) <- c('sample', 'gene', 'trv_type', 'label')
    } else {
        colnames(x) <- c('sample', 'gene', 'trv_type')
    }
    
    return(x)
}
#' Assign NA samples a gene
#'
#' Replace NA values in a gene column with the top gene name
#' @name waterfall_NA2gene
#' @param x a data frame in anno format
#' @return a data frame with NA values in a gene column coerced to the top gene
#' name
#' @noRd
#' @importFrom stats na.omit

waterfall_NA2gene <- function(x)
{
    # Get The gene with the most Mutations and add the NA samples to that gene
    # (ensures that the NAs are added in as gene with most mutations will always
    # be plotted) i.e. makes sure that samples are plotted, happens with rmvSilent param
    
    # find top gene
    top_gene <- stats::na.omit(rev(x$gene))[1]
    
    # set the trv_type to NA if gene is NA (makes sure that wile a sample is plotted the cell is empty)
    x$trv_type <- replace(x$trv_type, is.na(x$gene), NA)
    x$gene <- replace(x$gene, is.na(x$gene), top_gene)
    
    return(x)
}
#' align plots
#'
#' align mutation landscape, mutation burden on sample, and mutation burden on
#' gene plots
#' @name waterfall_align
#' @param heatmap ggplot object displaying a mutation landscape
#' @param genes ggplot object displaying mutation burden on gene
#' @param burden ggplot object displaying mutation burden on sample
#' @param clinical ggplot object displaying clinical information "optional"
#' @param proportion ggplot object displaying proportion of mutation types "optional"
#' @param section_heights Heights of each section (should sum to one)
#' @return a grob object
#' @noRd
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid nullGrob

waterfall_align <- function(genes, heatmap, burden, clinical, proportion, 
                            section_heights) {
    if (missing(section_heights)) {    
        if (missing(clinical)) {
            if (is.null(proportion)) {
                section_heights <- c(1, 4)
            } else {
                section_heights <- c(1, 4, 0.5)
            }
        } else {
            if (is.null(proportion)) {
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
    ## https://github.com/baptiste/gridextra/wiki/arranging-ggplot#legends
    ind_legend <- grep("guide", heatmap_grob$layout$name)
    heatmap_legend <- heatmap_grob[["grobs"]][[ind_legend]]
    heatmap_width <- sum(heatmap_legend$width)
    heatmap_grob <- ggplot2::ggplotGrob(heatmap + theme(legend.position="none"))
    
    if(grid::is.grob(burden)){
        burden_width <- NULL
        burden_grob <- burden
        blankPanel <- grid::grid.rect(gp=grid::gpar(col="white"))
        burden_legend <- grid::nullGrob()
    } else {
        burden_grob <- ggplot2::ggplotGrob(burden)
        ## Strip out legends and plot separately
        ind_legend <- grep("guide", burden_grob$layout$name)
        burden_legend <- burden_grob[["grobs"]][[ind_legend]]
        burden_width <- sum(burden_legend$width)
        burden_grob <- ggplot2::ggplotGrob(burden + theme(legend.position="none"))
        
        blankPanel <- grid::grid.rect(gp=grid::gpar(col="white"))
    }
    
    if (!is.null(proportion)) {
        prop_grob <- ggplot2::ggplotGrob(proportion)
        ind_legend <- grep("guide", prop_grob$layout$name)
        prop_legend <- prop_grob[["grobs"]][[ind_legend]]
        prop_width <- sum(prop_legend$width)
        prop_grob <- ggplot2::ggplotGrob(proportion + theme(legend.position="none"))
        
    } else prop_width <- NULL
    
    if(!missing(clinical)) {
        clin_grob <- ggplot2::ggplotGrob(clinical)
        ind_legend <- grep("guide", clin_grob$layout$name)
        clin_legend <- clin_grob[["grobs"]][[ind_legend]]
        clin_width <- sum(clin_legend$width)
        clin_grob <- ggplot2::ggplotGrob(clinical + theme(legend.position="none"))
        
    } else clin_width <- NULL
    
    ## https://github.com/baptiste/gridextra/wiki/arranging-ggplot#legends
    legend_width <- unit.pmax(heatmap_width, burden_width, clin_width, prop_width)
    width_left <- unit(1, "npc") - max(legend_width)
    widths <- grid::unit.c(width_left * 0.15, width_left * 0.85, legend_width)
    
    # Adjust the grob widths so heatmap and burden plots line up
    if(!missing(clinical)) {
        if (!is.null(proportion)) {
            maxwidth <- grid::unit.pmax(heatmap_grob$widths,
                                        burden_grob$widths,
                                        clin_grob$widths,
                                        prop_grob$widths)
            prop_grob$widths <- as.list(maxwidth)
        } else {
            maxwidth <- grid::unit.pmax(heatmap_grob$widths,
                                        burden_grob$widths,
                                        clin_grob$widths
            )
        }
        burden_grob$widths <- as.list(maxwidth)
        heatmap_grob$widths <- as.list(maxwidth)
        clin_grob$widths <- as.list(maxwidth)
    } else {
        if (!is.null(proportion)) {
            maxwidth <- grid::unit.pmax(heatmap_grob$widths,
                                        burden_grob$widths,
                                        prop_grob$widths)
            prop_grob$widths <- as.list(maxwidth)
        } else {   
            maxwidth <- grid::unit.pmax(heatmap_grob$widths,
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
        if (!is.null(proportion)) {
            nrows <- 4
            grobs <- list(
                blankPanel, burden_grob, burden_legend, 
                genes_grob, heatmap_grob, heatmap_legend, 
                blankPanel, prop_grob, prop_legend, 
                blankPanel, clin_grob, clin_legend)
        } else {
            nrows <- 3
            grobs <- list(
                blankPanel, burden_grob, burden_legend, 
                genes_grob, heatmap_grob, heatmap_legend,
                blankPanel, clin_grob, clin_legend)
        }
        heatmap <- gridExtra::arrangeGrob(grobs = grobs,
                                          ncol=3, nrow=nrows, 
                                          widths=widths,
                                          heights=section_heights)
        
        
        
    } else {
        if (!is.null(proportion)) {
            nrows <- 3
            grobs <- list(
                blankPanel, burden_grob, burden_legend,
                genes_grob, heatmap_grob, heatmap_legend, 
                blankPanel, prop_grob, prop_legend)
        } else {
            grobs <- list(
                blankPanel, burden_grob, burden_legend,
                genes_grob, heatmap_grob, heatmap_legend
            )
            nrows <- 2
        }
        heatmap <- gridExtra::arrangeGrob(grobs = grobs,
                                          ncol=3, nrow=nrows,
                                          widths=widths, heights=section_heights)
    }
    
    return(heatmap)
}
#' plot mutation recurrence in genes
#'
#' plot a bar graph displaying the percentage of samples with a mutation
#' @name waterfall_buildGenePrevelance
#' @param data_frame a data frame in MAF format
#' @param gene_label_size numeric value indicating the size of the gene labels
#'      on the y-axis
#' @param layers additional ggplot2 layers
#' @noRd
#' @return a ggplot object
#' @importFrom plyr count
#' @importFrom stats na.omit

waterfall_buildGenePrevelance <- function(data_frame, gene_label_size=8, layers=NULL)
{
    # Convert all silent mutations to Synonymous, and all else to non-synonymous
    data_frame$trv_type <- as.character(data_frame$trv_type)
    data_frame$trv_type[toupper(data_frame$trv_type) != toupper('silent')] <-
        'Non Synonymous'
    data_frame$trv_type[toupper(data_frame$trv_type) == toupper('silent')] <-
        'Synonymous'
    data_frame$trv_type <- factor(data_frame$trv_type,
                                  levels=c('Synonymous', 'Non Synonymous'))
    
    # Define the number of samples for the Percentage calculation
    # (Note: to pass a variable outside of aes into aes it needs to be
    # defined again)
    total_number_sample <- nlevels(data_frame$sample)
    
    # Make a count column
    data_frame <- plyr::count(data_frame, c("gene", "trv_type"))
    data_frame$prop <- data_frame$freq/total_number_sample * 100
    
    # Define Theme and various other layers to be passed to ggplot
    theme <- theme(axis.text.y=element_text(size=gene_label_size,
                                            colour='black', face='italic'),
                   axis.title.y=element_blank(),
                   legend.position=('none'))
    y_limits <- ylim(100, 0)
    y_label <- ylab('% Mutant')
    legend <- scale_fill_manual(name="Translational Effect",
                                values=c("Non Synonymous"="blue", "Synonymous"="red"),
                                breaks=c('Synonymous', 'Non Synonymous'),
                                drop=FALSE)
    
    if(!is.null(layers))
    {
        layers <- layers
    } else {
        layers <- geom_blank()
    }
    
    # Plotting call
    p1 <- ggplot(stats::na.omit(data_frame),
                 aes_string(x='gene', y='prop', fill='trv_type')) +
        geom_bar(position='stack', alpha=.75, width=1, stat='identity') +
        theme_bw() + coord_flip() + theme + y_label + scale_y_reverse() +
        legend + layers
    
    return(p1)
}
#' Plot a mutation heatmap
#'
#' Plot a Mutation Landscape with variables sample, gene, mutation
#' @name waterfall_buildMain
#' @param data_frame a data frame in MAF format
#' @param grid boolean value whether to overlay a grid on the plot
#' @param label_x boolean value whether to label the x axis
#' @param file_type character string specifying the file type, one of 'MAF' or
#' 'MGI'
#' @param drop_mutation Boolean specifying whether to drop unused
#' "mutation type" levels from the legend
#' @param plot_x_title Boolean specifying whether to plot the x_axis title
#' @param plot_label Boolean specifying whether to plot text inside each cell
#' @param plot_label_size Integer specifying text size of cell labels
#' @param plot_palette Character vector specifying colors to fill on mutation
#' type
#' @param layers additional ggplot2 layers to plot
#' @param plot_label_angle angle at which to plot label text if plot_label is
#' true
#' @noRd
#' @return a ggplot2 object
#' @import ggplot2

waterfall_buildMain <- function(data_frame, grid=TRUE, label_x=FALSE,
                                file_type='MGI',
                                drop_mutation=FALSE, plot_x_title=TRUE,
                                plot_label=FALSE, plot_label_size=4,
                                plot_palette=NULL, layers=NULL,
                                plot_label_angle=0)
{
    # NOTE: scale_x_discret(drop=FALSE), added to ggplot2 call to ensure samples
    # remove by 'mutation_recurrence_subset' are still plotted as empty tiles
    
    # Define layers
    
    # grid overlay
    vertical_grid <- geom_vline(xintercept = seq(.5, nlevels(data_frame$sample),
                                                 by=1),
                                linetype='solid', colour='grey80', size=.01)
    
    if(length(unique(data_frame$gene)) == 1)
    {
        horizontal_grid <- geom_blank()
    } else {
        horizontal_grid <- geom_hline(yintercept = seq(1.5, length(unique(data_frame$gene)),
                                                       by=1),
                                      linetype='solid', colour='grey80',
                                      size=.01)
    }
    
    # Declare the appropriate palette
    palette <- waterfall_select_palette(file_type = file_type, 
                                        custom_palette = plot_palette)
    breaks_labels <- waterfall_palette_names(palette,
                                             file_type, data_frame)
    breaks <- breaks_labels[["breaks"]]
    labels <- breaks_labels[["labels"]]
    
    if(drop_mutation == TRUE)
    {
        # Create Legend
        legend <- scale_fill_manual(name="Mutation Type", values=palette,
                                    breaks=breaks, labels=labels, drop=TRUE)
    } else if (drop_mutation == FALSE) {
        # Create Legend
        legend <- scale_fill_manual(name="Mutation Type", values=palette,
                                    breaks=breaks, labels=labels, drop=FALSE)
    }
    
    # X Label
    x_label <- xlab(paste0('Sample (n=', nlevels(data_frame$sample), ')'))
    
    if(plot_label == TRUE)
    {
        label <- geom_text(data=data_frame,
                           mapping=aes_string(x='sample', y='gene',
                                              label='label'),
                           size=plot_label_size, colour='white',
                           angle=plot_label_angle)
    } else {
        label <- geom_blank()
    }
    
    # Theme, Boolean, if specified to plot x labels, define theme such that
    # labels are plotted
    default_theme <- theme(axis.ticks=element_blank(),
                           panel.grid.major=element_blank(),
                           panel.grid.minor=element_blank(),
                           axis.title.y=element_blank(),
                           panel.background=element_rect(fill='white',
                                                         colour='white'),
                           axis.text.y=element_blank(),
                           plot.title=element_blank())
    if(label_x & plot_x_title)
    {
        theme <-  theme(axis.text.x=element_text(angle=50, hjust=1),
                        axis.title.x=element_text(size=10),
                        legend.title=element_text(size=14)
        )
    } else if(!label_x & plot_x_title) {
        theme <-  theme(axis.text.x=element_blank(),                        
                        axis.title.x=element_text(size=20),
                        legend.title=element_text(size=14),
                        panel.border=element_rect(colour='grey80', fill=NA,
                                                  size=.1),
                        legend.position=("right"))
    } else if(label_x & !plot_x_title) {
        theme <-  theme(axis.text.x=element_text(angle=50, hjust=1),
                        axis.title.x=element_blank(),
                        legend.title=element_text(size=14))
    } else if(!label_x & !plot_x_title) {
        theme <-  theme(
            axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            legend.title=element_text(size=14),
            panel.border=element_rect(colour='grey80', fill=NA,
                                      size=.1),
            legend.position=("right"))
    }
    
    # additional parameters
    if(!is.null(layers))
    {
        layers <- layers
    } else {
        layers <- geom_blank()
    }
    
    # ggplot call
    p1 <- ggplot(data_frame, aes_string('sample', 'gene')) +
        geom_tile(aes_string(fill='trv_type'), position="identity")
    if(grid == TRUE) {
        p1 <- p1 + vertical_grid + horizontal_grid 
    }
    p1 <- p1 +  default_theme + theme + legend + x_label +
        scale_x_discrete(drop=FALSE) + label + layers
    
    return(p1)
}
#' plot mutation burden
#'
#' plot a barchart showing mutations per MB
#' @name waterfall_buildMutBurden_A
#' @param x a data frame in MAF format
#' @param coverage_space an integer specifying the coverage space in base pairs
#' from which a mutation could occur
#' @param layers Additional ggplot2 layers to plot
#' @noRd
#' @return a ggplot object

waterfall_buildMutBurden_A <- function(x, coverage_space, layers=NULL)
{
    # Add in mutations per MB calculation
    x$mutation_per_MB <-
        x$mutation_total/coverage_space * 1000000
    
    # Alter GGplot2 Theme
    theme <- theme(axis.ticks.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   legend.title = element_text(size=14),
                   axis.text.y = element_text(colour = "black"),
                   axis.title.y = element_text(colour = "black"),
                   panel.background = element_blank(),
                   # panel.grid.minor.y = element_line(colour = "black"),
                   panel.grid.major.y = element_line(colour = "grey80"),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.x = element_blank(),
                   panel.border = element_rect(fill = NA)
    )
    
    # Add Legend
    legend <- scale_fill_manual(name="Translational Effect",
                                values=c("Synonymous"="red", "Non Synonymous"="blue"),
                                breaks=c("Synonymous", "Non Synonymous"),
                                drop=FALSE)
    
    # add y label
    y_label <- ylab('Mutations per MB')
    
    # ggplot2 call
    p1 <- ggplot(x, aes_string(x='sample', y='mutation_per_MB',
                               fill='trv_type')) +
        geom_bar(stat='identity', alpha=.75, width=1) +
        theme + y_label + legend + layers
    
    return(p1)
}
#' plot mutation burden
#'
#' plot a barchart showing mutation burden given by data frame
#' @name waterfall_buildMutBurden_B
#' @param x a data frame containing columns sample, mut_burden
#' @param layers additional ggplot2 layers to plot
#' @return a ggplot object
#' @noRd
#' @import ggplot2

waterfall_buildMutBurden_B <- function(x, layers=NULL)
{
    # add in fake column for legend
    # (necessary to have legend for proper plot alignment)
    x$Type <- c("Undefined")
    x$Type <- factor(x$Type,
                     levels=c("Synonymous", "Non Synonymous", "Undefined"))
    
    # Define Theme
    theme <- theme(axis.ticks.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   legend.title = element_text(size=14),
                   axis.text.y = element_text(colour = "black"),
                   axis.title.y = element_text(colour = "black"),
                   panel.background = element_blank(),
                   # panel.grid.minor.y = element_line(colour = "black"),
                   panel.grid.major.y = element_line(colour = "grey80"),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.x = element_blank(),
                   panel.border = element_rect(fill = NA)
    )
    
    # Define additional parameters
    y_label <- ylab("Mutation Burden")
    legend <- scale_fill_manual(name="Translational Effect",
                                values=c("Non Synonymous"="blue",
                                         "Synonymous"="red",
                                         "Undefined"="black"),
                                breaks=c("Synonymous", "Non Synonymous", "Undefined"),
                                drop=FALSE)
    
    
    bar <- geom_bar(stat='identity', alpha=.75, width=1)
    
    # ggplot2 call
    p1 <- ggplot(x, aes_string(x='sample', y='mut_burden', fill='Type')) + bar +
        theme + y_label + legend + layers
    
    return(p1)
}
#' @title Build mutational profile plot
#' 
#' @param data_frame input data.frame
#' @param plot_palette Color palette to use
#' @param file_type MAF, etc
#' @param layers layer(s) to add to this plot object
#' @param x_label label this plot?
#' @description Builds a ggplot object showing individuals' mutational profile
#' @noRd
#' @return a ggplot object
waterfall_build_proportions <- function(data_frame, plot_palette, 
                                        file_type, layers, x_label) {
    
    # Declare the appropriate palette
    palette <- waterfall_select_palette(file_type, custom_palette = plot_palette)
    breaks_labels <- waterfall_palette_names(palette, file_type, data_frame)
    breaks <- breaks_labels[["breaks"]]
    labels <- breaks_labels[["labels"]]
    
    if (x_label) {
        x_label_obj <- xlab(paste0('Sample (n=', nlevels(data_frame$sample), ')'))
    } else {
        x_label_obj <- geom_blank()
    }
    
    p5 <- ggplot(data_frame, aes_string(x = 'sample', fill = 'trv_type')) + 
        geom_bar(position = "fill", width = 0.95) + 
        scale_y_continuous(
            name = "Proportion", 
            labels = scales::percent_format()) + 
        scale_fill_manual(name="Mutation Type", 
                          values=palette, 
                          breaks = breaks,
                          labels = labels
        ) + 
        guides(fill = guide_legend(title = "Mutation Type", ncol = 2)) +
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x =  element_blank(),
            axis.title.x = if (x_label) element_text() else element_blank(),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(colour = "black", hjust = 0),
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank()
        ) + layers + x_label_obj
    p5
}
#' Calculate Synonymous/Nonsynonymous mutation frequency
#'
#' Creates a data frame giving synonymous/nonsynonymous counts on a sample level
#' @name waterfall_calcMutFreq
#' @param x data frame in long format with columns sample, trv_type
#' @importFrom data.table melt
#' @noRd
#' @return a data frame with synonymous/nonsynonymous counts appended

waterfall_calcMutFreq <- function(x)
{
    message("Calculating frequency of mutations...")
    # Change trv_type calls to either synonymous or non synonymous,
    # for use in the mutation per Mb plot
    x$trv_type <- as.character(x$trv_type)
    x$trv_type[toupper(x$trv_type) != toupper('silent')] <- 'Non Synonymous'
    x$trv_type[toupper(x$trv_type) == toupper('silent')] <- 'Synonymous'
    x$trv_type <- factor(x$trv_type, levels=c('Synonymous', 'Non Synonymous'))
    
    # Obtain a data frame of mutation counts on the sample level
    mutation_counts <- table(x[,c('sample', 'trv_type')])
    mutation_counts <- as.data.frame(data.table::melt(mutation_counts))
    colnames(mutation_counts) <- c('sample', 'trv_type', 'mutation_total')
    
    return(mutation_counts)
}
#' mutation sample cutoff gene based
#'
#' Subset a internal mutSpec file keeping only samples within the specified gene
#'  list
#' @name waterfall_geneAlt
#' @param x a data frame in long format with columns 'gene', 'trv_type'
#' @param genes character vector listing genes to plot
#' @noRd
#' @return a subset data frame

waterfall_geneAlt <- function(x, genes)
{
    message("Removing genes not in: ", toString(genes))
    # Perform quality checks
    if(typeof(genes) != 'character' & class(genes) != 'character')
    {
        memo <- paste0("argument supplied to plotGenes is not a character ",
                       "vector, attempting to coerce")
        warning(memo)
        genes <- as.character(genes)
    }
    
    if(!all(toupper(genes) %in% toupper(x$gene)))
    {
        memo <- paste0("genes supplied in plotGenes contains an element not ",
                       "found in x or it's subsequent subsets")
        warning(memo)
    }
    genes <- c(genes, NA)
    x <- x[(toupper(x$gene) %in% toupper(genes)), ]
    
    return(x)
}
#' Mutation Recurrence Cutoff
#'
#' Subset a MAF file keeping only samples that meet a mutation recurrence cutoff
#' @name waterfall_geneRecurCutoff
#' @param x data frame in long format with columns 'gene', 'trv_type', 'sample'
#' @param recurrence_cutoff integer specifying removal of entries not seen
#' in at least "x" percent of samples
#' @return a subset data frame
#' @noRd
#' @importFrom plyr count
#' @importFrom stats na.omit

waterfall_geneRecurCutoff <- function(x, recurrence_cutoff)
{
    if(recurrence_cutoff != 0)
    {
        message("Performing recurrence cutoff...")
    }
    mutRecur <- plyr::count(unique(x), vars=c("gene"))
    mutRecur <- stats::na.omit(mutRecur)
    mutRecur$prop <- mutRecur$freq/nlevels(x$sample)
    
    # If recurrence cutoff specified exceeds upper limit such that no plot
    # useful would be generated, reset recurrence cutoff
    maxRecur <- max(mutRecur$prop)
    if(maxRecur < recurrence_cutoff)
    {
        memo <- paste0("The recurrence cutoff specified exceeds the recurrence",
                       " seen in the data, resetting this value to equal max ",
                       "recurrence:", maxRecur)
        warning(memo)
        recurrence_cutoff <- maxRecur
    }
    
    gene_above_recur <- mutRecur[mutRecur$prop >= recurrence_cutoff,]$gene
    # add NA to the end of 'gene_above_recurrence' vector, allowing for all
    # samples having NA as a gene name to be retained in the subset below
    gene_above_recur <- c(as.character(gene_above_recur), NA)
    
    # subset the original data frame based on the following: keep gene if it is
    # in the gene vector in "mutation_recurrence_subset"
    x <- x[(x$gene %in% gene_above_recur), ]
    
    return(x)
}
#' sort waterfall file by gene
#'
#' order a waterfall file ranking genes with more mutations higher if a gene
#' order is unspecified.
#' @name waterfall_geneSort
#' @param x Data frame with columns names "gene", "trv_type".
#' @param geneOrder Character vector specifying the order in which to plot
#' genes.
#' @noRd
#' @return Character vector of ordered genes

waterfall_geneSort <- function(x, geneOrder=NULL)
{
    if(!is.null(geneOrder))
    {
        geneOrder <- as.character(unique(geneOrder))
        # if there are any genes in geneOrder not in x, remove those
        gene_to_rmv <- geneOrder[!geneOrder %in% unique(x$gene)]
        if(length(gene_to_rmv) == 0 & length(geneOrder) >= 1)
        {
            return(rev(geneOrder))
        } else if(length(gene_to_rmv) == length(geneOrder)) {
            memo <- paste0("Did not find any genes supplied to the parameter",
                           " geneOrder in the input supplied to x, perhaps ",
                           " there is a case issue?")
            warning(memo)
        } else if(length(gene_to_rmv) != length(geneOrder)) {
            memo <- paste0("The following genes were not found in the input", 
                           " supplied to parameter x: ", toString(gene_to_rmv),
                           ", removing these from geneOrder!")
            warning(memo)
            gene_order <- geneOrder[geneOrder %in% unique(x$gene)]
            return(rev(gene_order))
        }
    }
    
    # order based on mutation frequency
    gene_mutation_table <- table(x[,c('gene', 'trv_type')])
    gene_order <- names(sort(rowSums(gene_mutation_table)))
    return(gene_order)
}
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
#' @noRd
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
        } else if(toupper(file_type) == toupper('Custom')) {
            memo <- paste0("Detected NULL in variant_class_order, ",
                           "this parameter is required if file_type is set ",
                           "to \"Custom\"")
            stop(memo)
        }
    } else {
        mutation_order <- unique(c(variant_class_order, NA))
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
#'  waterfall_palette_names
#' 
#' @title waterfall_palette_names
#' 
#' @description Make labels and breaks for palettes
#' 
#' @param palette Named colour vector as input 
#' @param file_type Which file type is involved?
#' @param data_frame Only used if file_type is "custom"
#' @noRd
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

#' Check input to mutSpec
#'
#' Perform a data quality check on input to mutSpec
#' @name waterfall_qual
#' @param x a data frame in annotation format
#' @param y a data frame containing clinical data or a null object
#' @param z a data frame containing mutation burden information or a null object
#' @param file_type Character string specifying the input format to expect in x
#' @param label_col Character string specifying the column name of a label
#' column
#' @noRd
#' @return a list of data frames passing quality checks

waterfall_qual <- function(x, y, z, file_type, label_col)
{
    # print message statement
    message("Checking if input is properly formatted...")
    
    # Check input data to x
    if(!is.data.frame(x))
    {
        stop("Did not detect a data frame for input to x")
    }
    
    # Convert file type to internal format
    if(toupper(file_type) == toupper("MAF"))
    {
        x <- waterfall_MAF2anno(x, label_col)
    } else if(toupper(file_type) == toupper("MGI")) {
        x <- waterfall_MGI2anno(x, label_col)
    } else if(toupper(file_type) == toupper("Custom")) {
        x <- waterfall_Custom2anno(x, label_col)
    } else {
        stop("Unrecognized file_type: ", file_type)
    }
    
    # drop unused levels in x
    x$sample <- as.factor(x$sample)
    x$gene <- as.factor(x$gene)
    x$trv_type <- as.factor(x$trv_type)
    x <- droplevels(x)
    
    # Check input data to clinDat
    if(!is.null(y))
    {
        if(!is.data.frame(y))
        {
            stop("Did not detect a data frame for input to clinDat")
        }
        y <- droplevels(y)
        
        if(!all(c('sample', 'variable', 'value') %in% colnames(y)))
        {
            stop("Did not detect correct sample names in clinDat")
        }
        
        # make sure clinical data columns are of expected class
        y$sample <- as.factor(y$sample)
        y$variable <- as.factor(y$variable)
        y$value <- as.character(y$value)
    }
    
    # check input data to mutBurden
    if(!is.null(z))
    {
        if(!is.data.frame(z))
        {
            stop("Did not detect a data frame for input to mutBurden")
        }
        z <- droplevels(z)
        
        if(!all(c('sample', 'mut_burden') %in% colnames(z)))
        {
            stop("Did not detect correct sample names in mutBurden")
        }
        
        # Make sure mutation burden columns are of proper class
        z$sample <- as.factor(z$sample)
        z$mut_burden <- as.numeric(as.character(z$mut_burden))
    }
    
    return(list(x, y, z))
}
#' Silent Mutation Removal
#'
#' Subset a MAF file setting keeping only sample information if a mutation
#' is silent
#' @name waterfall_rmvSilent
#' @param x a data frame with columns 'sample', 'gene', 'trv_type'
#' @noRd
#' @return a subset data frame

waterfall_rmvSilent <- function(x)
{
    message("Removing silent mutations...")
    
    # Index and remove those rows which contain silent mutations
    x[which(toupper(x$trv_type) == toupper('silent')), c('gene')] <- NA
    if("label" %in% colnames(x)){
        x[which(toupper(x$trv_type) == toupper('silent')), c('label')] <- NA
    }  
    x[which(toupper(x$trv_type) == toupper('silent')), c('trv_type')] <- NA
    
    return(x)
}
#' mutation sample subset sample based
#'
#' Alter a mutSpec input file keeping/adding entries in a selection of samples
#' @name waterfall_sampAlt
#' @param x a data frame in long format with columns 'sample', 'trv_type'
#' @param samples character vector giving samples to plot
#' @return a subset data frame
#' @noRd
#' @importFrom plyr rbind.fill

waterfall_sampAlt <- function(x, samples)
{
    message("Retrieving requested samples from supplied data...")
    # Perform quality check
    if(typeof(samples) != 'character' & class(samples) != 'character')
    {
        memo <- paste0("argument supplied to main.samples is not a ",
                       "character vector, attempting to coerce")
        warning(memo)
        samples <- as.character(samples)
    }
    
    # add in samples not in the original data frame x
    if(!all(toupper(samples) %in% toupper(x$sample)))
    {
        memo <- paste0("Detected one or more samples supplied to main.samples ",
                       "not found in x ... adding samples to plot")
        message(memo)
        
        samp_not_in_x <- as.data.frame(samples[which(!(samples %in% x$sample))])
        colnames(samp_not_in_x) <- "sample"
        x <- plyr::rbind.fill(samp_not_in_x, x)
    }
    
    # Subset the data based on the arguments supplied in main.samples,
    # then relevel the data frame
    x <- x[(toupper(x$sample) %in% toupper(samples)), ]
    x <- droplevels(x)
    x$sample <- factor(x$sample, levels=samples)
    
    return(x)
}
#' sort samples in an internal waterfall file.
#'
#' perform a hiearchial sort on samples based on the presence of mutations in
#' an ordered list of genes if a sample order is unspecified.
#' @name waterfall_sampSort
#' @param x a data frame in long format with column names "sample",
#' "gene", "trv_type"
#' @param sampOrder Character vector specifying the order of samples to plot.
#' @return a vector of samples in a sorted order
#' @noRd
#' @importFrom data.table dcast

waterfall_sampSort <- function(x, sampOrder=NULL)
{
    # if a sample order is already defined plot that instead
    if(!is.null(sampOrder))
    {
        sampOrder <- as.character(unique(sampOrder))
        
        # determine if there are any new samples in sampOrder not in the data
        # and remove if true
        new_samples <- sampOrder[!sampOrder %in% levels(x$sample)]
        if(length(new_samples) != 0)
        {
            memo <- paste0("The following samples were not detected in the ",
                           "original data: ", toString(new_samples),
                           ", removing these samples from sampOrder!")
            warning(memo)
            sampOrder <- sampOrder[!sampOrder %in% new_samples]
        }
        
        # return what was given originally
        return(sampOrder)
    }
    
    # recast the data going from long format to wide format, values in this data
    # are counts of a mutation call
    wide_data <- data.table::dcast(x, sample ~ gene, fun.aggregate = length,
                                 value.var="trv_type")
    
    # apply a boolean function to convert the data frame values to 1's and 0's
    values <- wide_data[,-1, drop=FALSE]
    sample <- wide_data[,1]
    values <- data.frame(apply(values, 2,
                               function(x) as.numeric(as.logical(x))))
    wide_boolean <- cbind(sample, values)
    
    # reverse the columns so that genes with highest mutation's are listed first
    # (assumes gene_sort has been run on the data frame)
    wide_boolean <- wide_boolean[,c(1, rev(2:ncol(wide_boolean)))]
    
    # if there are any NA values present in a sample at the gene put that
    # remove that sample and save (at this stage it should be samples)
    if(any(grepl("^NA.$", colnames(wide_boolean))))
    {
        # Find which column has the NA header
        NA_index <- which(grepl("^NA.$", colnames(wide_boolean)))
        
        # Append NA column to end of data frame
        NA_gene <- wide_boolean[,NA_index]
        wide_boolean <- wide_boolean[,-NA_index]
        
        # Save copy and remove samples with no mutations,
        # these will be added to the end
        samp_no_mut <- wide_boolean[rowSums(wide_boolean[2:ncol(wide_boolean)]) == 0,]$sample
        samp_no_mut <- as.character(samp_no_mut)
        wide_boolean <- wide_boolean[!wide_boolean$sample %in% samp_no_mut,]
    } else {
        samp_no_mut <- NULL
    }
    
    # hiearchial sort on all column's (i.e. genes) such that samples are
    # rearranged if there is a mutation in that gene
    sample_order <- wide_boolean[do.call(order, as.list(-wide_boolean[2:ncol(wide_boolean)])),]$sample
    
    # Put those samples not in sample order in from the original levels of the
    # data (these are samples with no mutations)
    not_in <- as.character(levels(sample_order)[which(!levels(sample_order) %in% sample_order)])
    not_in <- not_in[!not_in %in% samp_no_mut]
    sample_order <- c(as.character(sample_order), as.character(not_in))
    
    # Put those samples with no mutations back in
    if(!is.null(samp_no_mut))
    {
        sample_order <- c(sample_order, samp_no_mut)
    }
    
    return(sample_order)
}
#' waterfall_select_palette
#' 
#' @title Helper function to select a colour palette
#' 
#' @param file_type Which file tyoe is used?
#' @param custom_palette Nullable custom colour palette
#' @noRd
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

#' Construct transition-transversion plot
#'
#' Given a data frame construct a plot displaying the proportion or frequency of
#' transition and transversion types observed in a cohort.
#' @name TvTi
#' @param x Object of class data frame with rows representing transitions and
#' transversions. The data frame must contain the following columns 'sample',
#' reference' and 'variant' or alternatively "Tumor_Sample_Barcode",
#' "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2" depending on the
#' argument supplied to the fileType parameter. (required)
#' @param y Named vector or data frame representing the expected transition and 
#' transversion rates. Either option must name transition and transverions as
#' follows: "A->C or T->G (TV)", "A->G or T->C (TI)", "A->T or T->A (TV)",
#'  "G->A or C->T (TI)", "G->C or C->G (TV)", "G->T or C->A (TV)". If specifying
#'  a data frame, the data frame must contain the following columns names
#'  "Prop", "trans_tranv" (optional see vignette).
#' @param clinData Object of class data frame with rows representing clinical
#' data. The data frame should be in "long format" and columns must be names as
#' "sample", "variable", and "value" (optional see details and vignette).
#' @param type Character string specifying if the plot should display the
#' Proportion or Frequency of transitions/transversions observed. One of
#' "Proportion" or "Frequency", defaults to "Proportion".
#' @param lab_Xaxis Boolean specifying whether to label the x-axis in the plot.
#' @param lab_txtAngle Integer specifying the angle of labels on the x-axis of
#' the plot.
#' @param palette Character vector of length 6 specifying colours for each
#' of the six possible transition transversion types.
#' @param fileType Character string specifying the format the input given to
#' parameter x is in, one of 'MAF', 'MGI'. The former option requires the data
#' frame given to x to contain the following column names
#' "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele1",
#' "Tumor_Seq_Allele2" the later option requires the data frame givin to x to
#' contain the following column names "reference", "variant" and "sample".
#' (required)
#' @param tvtiLayer Valid ggplot2 layer to be added to the main plot.
#' @param expecLayer Valid ggplot2 layer to be added to the expected sub-plot.
#' @param sort Character string specifying the sort order of the sample
#' variables in the plot. Arguments to this parameter should be "sample",
#' "tvti", or "none" to sort the x-axis by sample name, transition transversion
#' frequency, or no sort respectively.
#' @param out Character vector specifying the the object to output, one of
#' "data", "grob", or "plot", defaults to "plot" (see returns).
#' @param clinLegCol Integer specifying the number of columns in the legend for
#' the clinical data, only valid if argument is supplied to parameter clinData.
#' @param clinVarCol Named character vector specifying the mapping of colours
#' to variables in the variable column of the data frame supplied to clinData
#' (ex. "variable"="colour").
#' @param clinVarOrder Character vector specifying the order in which to plot
#' variables in the variable column of the argument given to the parameter
#' clinData. The argument supplied to this parameter should have the same unique
#' length and values as in the variable column of the argument supplied to 
#' parameter clinData (see vignette).
#' @param clinLayer Valid ggplot2 layer to be added to the clinical sub-plot.
#' @param progress Boolean specifying if progress bar should be displayed for
#' the function.
#' @param sample_order_input Sample orders to be used
#' @param layers ggplot object to be added to proportions plot
#' @param return_plot Return as ggplot object? Only returns main plot
#' @details TvTi is a function designed to display proportion or frequency
#' of transitions and transversion seen in a data frame supplied to parameter x.
#' @examples
#' TvTi(brcaMAF, type='Frequency',
#' palette=c("#77C55D", "#A461B4", "#C1524B", "#93B5BB", "#4F433F", "#BFA753"),
#' lab_txtAngle=60, fileType="MAF")
#' @return One of the following, a list of dataframes containing data to be
#' plotted, a grob object, or a plot.
#' @importFrom plyr adply
#' @importFrom gtools mixedsort
#' @export

TvTi <- function(x, fileType=NULL, y=NULL, clinData=NULL, type='Proportion',
                 lab_Xaxis=TRUE, lab_txtAngle=45,
                 palette=c('#D53E4F', '#FC8D59', '#FEE08B', '#E6F598',
                           '#99D594', '#3288BD'),
                 tvtiLayer=NULL, expecLayer=NULL,
                 sort='none', clinLegCol=NULL, clinVarCol=NULL,
                 clinVarOrder=NULL, clinLayer=NULL, progress=TRUE, out="plot",
                 sample_order_input, layers = NULL, return_plot = FALSE)
{ 
    message("This function has been deprecated in order to implement an object oriented programming style! Please use MutSpectra() instead!")    
    # Perform quality checks
    output <- TvTi_qual(x, y, clinData, file_type=fileType)
    x <- output$input1
    y <- output$input2
    clinData <- output$input3
    
    # add transition/transversion info
    if(isTRUE(progress))
    {
        message("annotating transitions and transversions")
        x <- plyr::adply(x, 1, TvTi_annoTransTranv, .progress='text')
    } else {
        x <- plyr::adply(x, 1, TvTi_annoTransTranv)
    }
    
    
    # Calculate the proportion of transitions/transversions
    x <- TvTi_calcTransTranvFreq(x)
    
    # re-level based on proportion values or via a smart sort or not at all
    if(toupper(sort) == toupper('sample'))
    {
        sample_order <- as.vector(unique(x$sample))
        sample_order <- gtools::mixedsort(sample_order)
        x$sample <- factor(x$sample, levels=sample_order)
    } else if(toupper(sort) == toupper('tvti')) {
        sample_order <- x[order(x$trans_tranv, -x$Prop),]
        sample_order <- sample_order[sample_order$Prop != 0,]
        sample_order <- unique(sample_order$sample)
        x$sample <- factor(x$sample, levels=sample_order)
    } else if(toupper(sort) == toupper("custom")) {
        sample_order <- sample_order_input
        x$sample <- factor(x$sample, levels = sample_order)
    } else if(toupper(sort) == toupper('none')){
        sample_order <- levels(x$sample)
    } else {
        memo <-paste0(sort, " is not a valid parameter for sort, please",
                      " specify one of \"sample\", \"tvti\", \"custom\", \"none\"")
        stop(memo)
    }
    
    # Perform a quality control on y to ensure fill levels match x
    if(!is.null(y))
    {
        y$trans_tranv <- factor(y$trans_tranv, levels=levels(x$trans_tranv))
    }
    
    # Build the Transition/Transversion Plot if clinical data does not exist
    if(is.null(clinData))
    {
        p1 <- TvTi_buildMain(x, y, type=type,
                             x_axis_text_angle=lab_txtAngle,
                             palette=palette, label_x_axis=lab_Xaxis,
                             tvti.layers=tvtiLayer, expec.layers=NULL,
                             title_x_axis = TRUE)
    }
    
    
    # Plot Clinical Data if Speccified and build modified TvTi main plot
    if(!is.null(clinData))
    {
        clinData$sample <- factor(clinData$sample, levels=sample_order)
        p3 <- multi_buildClin(clinData, clin.legend.col=clinLegCol, 
                              clin.var.colour=clinVarCol, 
                              clin.var.order=clinVarOrder,
                              clin.layers=clinLayer)
        
        # Build the Transition/Transversion Plot
        p1 <- TvTi_buildMain(x, y, type=type,
                             x_axis_text_angle=lab_txtAngle,
                             palette=palette, label_x_axis=lab_Xaxis,
                             tvti.layers=tvtiLayer, expec.layers=NULL,
                             title_x_axis=FALSE)
    } else {
        p3 <- NULL
    }
    
    if(!is.null(y))
    {
        # If y is input plot the expected values
        p2 <- TvTi_buildMain(y, y, type=type,
                             x_axis_text_angle=lab_txtAngle,
                             palette=palette, label_x_axis=lab_Xaxis,
                             plot_expected=TRUE, tvti.layers=NULL,
                             expec.layers=expecLayer)
    } else {
        p2 <- NULL
    }
    p1 <- p1 + layers
    # Decide what to output
    if (return_plot) {
        return(p1)
    } else {
        finalPlot <- TvTi_alignPlot(p1=p1, p2=p2, p3=p3)
        dataOut <- list("main"=x, "expect"=y)
        output <- multi_selectOut(data=dataOut, plot=finalPlot, draw=TRUE, out=out)
        return(output)
    }
}

#' align TvTi plots on y axis
#'
#' align transition/transversion plots
#' @name TvTi_alignPlot
#' @param p1 main plot
#' @param p2 left expected value subplot
#' @param p3 bottom clinical subplot
#' @noRd
#' @return ggplot object

TvTi_alignPlot <- function(p1=NULL, p2=NULL, p3=NULL)
{
    # define the ggplot's as grobs and create a blank plot
    gA <- ggplot2::ggplotGrob(p1)
    
    # Adjust the grob heights so p1, and p2 plots line up if p2 exists
    if(!is.null(p2))
    {
        # convert expected plot to grob
        gB <- ggplot2::ggplotGrob(p2)
        
        maxheight = grid::unit.pmax(gA$heights, gB$heights)
        gA$heights <- as.list(maxheight)
        gB$heights <- as.list(maxheight)
    }
    
    # adjust the grob widths so p1 and p3 line up if p3 exists
    if(!is.null(p3))
    {
        # Convert clinical plot to grob
        gC <- ggplot2::ggplotGrob(p3)
        gD <- grid::grid.rect(gp=grid::gpar(col="white"))
        
        maxwidth = grid::unit.pmax(gA$widths, gC$widths)
        gA$widths <- as.list(maxwidth)
        gC$widths <- as.list(maxwidth)        
    }
    
    # Build the final plot
    if(is.null(p3) & is.null(p2))
    {
        finalPlot <- gridExtra::arrangeGrob(gA, ncol=1, nrow=1)
    } else if(is.null(p3) & !is.null(p2)) {
        finalPlot <- gridExtra::arrangeGrob(gB, gA, ncol=2, nrow=1, widths=c(1,6))
    } else if(!is.null(p3) & is.null(p2)) {
        finalPlot <- gridExtra::arrangeGrob(gA, gC, ncol=1, nrow=2, heights=c(6,1))
    } else if(!is.null(p3) & !is.null(p2)) {
        finalPlot <- gridExtra::arrangeGrob(gB, gA, gD, gC, ncol=2, nrow=2, heights=c(6,1), widths=c(1,6))
    }
    
    return(finalPlot)
}

#' Annotate Transitions and Transversions
#'
#' Given a data frame with columns reference and variant annotate the base
#' change occurring
#' @name TvTi_annoTransTranv
#' @param x Object of class data frame containing columns 'reference', 'variant'
#' @noRd
#' @return Object of class data frame with transition/transversion annotations
#' appended

TvTi_annoTransTranv <- function(x)
{
    # add an extra column with the reference to variant
    x$base_change <- paste0(toupper(x$reference), "2", toupper(x$variant))
    
    # annotate the grouping of the base change
    x$trans_tranv <- switch(x$base_change, A2C="A->C or T->G (TV)",
                            T2G="A->C or T->G (TV)", A2G="A->G or T->C (TI)",
                            T2C="A->G or T->C (TI)", A2T="A->T or T->A (TV)",
                            T2A="A->T or T->A (TV)", G2A="G->A or C->T (TI)",
                            C2T="G->A or C->T (TI)", G2C="G->C or C->G (TV)",
                            C2G="G->C or C->G (TV)", G2T="G->T or C->A (TV)",
                            C2A="G->T or C->A (TV)")
    
    # remove the temp base change column
    x$base_change <- NULL
    return(x)
}

#' build transitions/transversions
#'
#' Given a data frame with columns 'trans_tranv', 'sample', 'Freq', and 'Prop',
#' build a transition/transversion plot
#' @name TvTi_buildMain
#' @param x Object of class data frame containing columns 'trans_tranv',
#' 'sample', 'Freq', and 'Prop'
#' @param y Object of class data frame containing columns 'Prop', 'trans_tranv'
#' for display of expected results
#' @param type Object of class character specifying whether to plot the
#' Proportion or Frequency, one of "Prop"
#' @param label_x_axis boolean specifying wheter to label x axis
#' @param x_axis_text_angle Integer specifying the angle to labels on x_axis
#' @param palette Character vector of length 6 specifying colors for
#' trans/tranv type
#' @param plot_expected Boolean specifying if this is the main TvTi plot or a
#' sub plot for expected values
#' @param tvti.layers Additional ggplot2 layers for the main plot
#' @param expec.layers Additional ggplot2 layers for the expected values plot
#' @param title_x_axis boolean specifying whether to display an x axis title
#' @noRd
#' @return GGplot Object
#' @import ggplot2

TvTi_buildMain <- function(x, y=NULL, type='Proportion', label_x_axis=TRUE,
                           x_axis_text_angle=45,
                           palette=c('#D53E4F', '#FC8D59', '#FEE08B', '#E6F598',
                                     '#99D594', '#3288BD'),
                           plot_expected=FALSE, tvti.layers=NULL,
                           expec.layers=NULL, title_x_axis=TRUE)
{
    
    if(!is.null(y))
    {
        # cumulativley sum the expected values and plot
        y$cumsum <- cumsum(y$Prop)
        expected <- geom_hline(data=y, mapping=aes_string(yintercept='cumsum'),
                               linetype="longdash", size=.5)
    } else {
        expected <- geom_blank()
    }
    
    # Define various parameters of plot
    if(plot_expected == TRUE)
    {
        bar <- geom_bar(data=x, mapping=aes_string(x=shQuote('Expected'),
                                                   y='Prop',
                                                   fill='trans_tranv'),
                        stat='identity', width=1)
    } else if(toupper(type) == 'PROPORTION') {
        bar <- geom_bar(data=x,
                        mapping=aes_string(x='sample', y='Prop',
                                           fill='trans_tranv'),
                        stat='identity', width=1)
    } else if(toupper(type) == 'FREQUENCY') {
        bar <- geom_bar(data=x,
                        mapping=aes_string(x='sample', y='Freq',
                                           fill='trans_tranv'),
                        stat='identity', width=1)
    }
    
    ylabel <- ylab(type)
    
    if(title_x_axis == TRUE)
    {
        xlabel <- xlab(paste0("Sample: n=", length(unique(x$sample))))
    } else {
        xlabel <- xlab('')
    }
    
    fill_palette <- scale_fill_manual(name='Transition/Transversion',
                                      values=palette)
    
    # additional layers to plot?
    if(!is.null(tvti.layers))
    {
        layers <- tvti.layers
    } else {
        layers <- geom_blank()
    }
    
    # Define theme of plot
    if(plot_expected == TRUE)
    {
        theme <- theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       legend.position='none',
                       axis.title.x=element_blank(),
                       axis.text.x=element_text(angle=x_axis_text_angle,
                                                hjust=1, vjust=1))
        if(!is.null(expec.layers))
        {
            layers <- expec.layers
        } else {
            layers <- geom_blank()
        }
        
    } else if(label_x_axis == TRUE) {
        theme <- theme(axis.text.x=element_text(angle=x_axis_text_angle,
                                                hjust=1, vjust=1))
    } else {
        theme <- theme(axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
    }
    
    # Define plot
    tmp <- data.frame(x=0, y=0)
    p1 <- ggplot(data=tmp, aes(y=0)) + bar + xlabel + ylabel +
        theme_bw() + theme + fill_palette + expected +
        guides(fill=guide_legend(reverse=TRUE)) + layers
    
    return(p1)
}

#' Calculate Transition/Transversion Frequency
#'
#' Given a data frame with columns reference, variant, sample, and trans/tranv
#' calculate the frequencies of transitions and transversion occuring.
#' @name TvTi_calcTransTranvFreq
#' @param x Object of class data frame containing columns 'reference',
#' 'variant', 'sample', 'trans_tranv'
#' @noRd
#' @return Object of class data frame with Frequency and Proportion of
#' Transistions/Transversions appended on a sample level

TvTi_calcTransTranvFreq <- function(x)
{
    # Ensure all possible combinations of trans/tranv are represented
    trans_tranv <- c("A->C or T->G (TV)", "A->G or T->C (TI)",
                     "A->T or T->A (TV)", "G->A or C->T (TI)",
                     "G->C or C->G (TV)", "G->T or C->A (TV)")
    sample <- c('dummy_sample')
    reference <- c('A')
    variant <- c('T')
    dummy_data <- data.frame(reference, variant, sample, trans_tranv)
    x <- rbind(dummy_data, x)
    
    # calculate the frequency of transitions/transversions on a sample basis
    x_freq <-  table(x$trans_tranv, x$sample)
    
    # calculate the proportion of transitions/transversions on a sample basis
    x_prop <-  prop.table(x_freq, 2)
    
    # format and remove the dummy data introduced above
    x_freq <- as.data.frame(x_freq)
    x_prop <- as.data.frame(x_prop)
    x <- cbind(x_freq, x_prop$Freq)
    colnames(x) <- c('trans_tranv', 'sample', 'Freq', 'Prop')
    x <- x[which(x$sample != "dummy_sample"),]
    
    return(x)
}

#' Convert .maf format to internal format
#'
#' Convert data frame in .maf format to an internally recogized format
#' @name TvTi_convMAF
#' @param x Object of class data frame containing columns
#' 'Tumor_Sample_Barcode', 'Reference_Allele' 'Tumor_Seq_Allele1',
#' 'Tumor_Seq_Allele2'
#' @noRd
#' @return a data frame, with column names 'sample', 'reference', 'variant'

TvTi_convMaf <- function(x)
{
    # Take out the appropriate columns and format for each allele
    x <- x[,c('Tumor_Sample_Barcode', 'Reference_Allele', 'Tumor_Seq_Allele1',
              'Tumor_Seq_Allele2')]
    
    allele1 <- x[,c('Tumor_Sample_Barcode', 'Reference_Allele',
                    'Tumor_Seq_Allele1')]
    colnames(allele1) <- c('sample', 'reference', 'variant')
    allele2 <- x[,c('Tumor_Sample_Barcode', 'Reference_Allele',
                    'Tumor_Seq_Allele2')]
    colnames(allele2) <- c('sample', 'reference', 'variant')
    
    #!!!Developer note: The if's are here because subsetting when there is
    #!!! nothing to subset (integer0) causes problems
    
    # if the tumor allele 1 matchest tumor allele2 remove that information from
    # one of the alleles
    if(any(as.character(allele1$variant) == as.character(allele2$variant)))
    {
        allele1 <- allele1[-which(as.character(allele1$variant) == as.character(allele2$variant)),]
    }
    
    # if the allele matches the reference remove it from the data
    if(any(as.character(allele1$reference) == as.character(allele1$variant)))
    {
        allele1 <- allele1[-which(as.character(allele1$reference) == as.character(allele1$variant)),]
    }
    if(any(as.character(allele2$reference) == as.character(allele2$variant)))
    {
        allele2 <- allele2[-which(as.character(allele2$reference) == as.character(allele2$variant)),]
    }
    
    # bind the data from both alleles together
    x <- rbind(allele1, allele2)
    
    return(x)
}

#' Check input to TvTi
#'
#' Perform quality check for input to function TvTi
#' @name TvTi_qual
#' @param x Object of class data frame containing columns 'sample', reference',
#' 'variant' for 'MGI' file or 'Tumor_Sample_Barcode', 'Reference_Allele',
#' 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2' for 'MAF' file
#' @param y Object of class data frame containing columns "Prop", "trans_tranv"
#' @param z Object of class data frame containing columns "sample", "variable",
#' "value" denoting clinical information
#' @param file_type Character string spedifying th input file type expected
#' @noRd
#' @return a data frame, or list of data frames passing quality checks

TvTi_qual <- function(x, y=NULL, z=NULL, file_type='MAF')
{
    # Check file type is valid
    if(!grepl("MAF|MGI", file_type))
    {
        memo <- paste0("Did not recognize input to paramter fileType as a",
                       " valid argument... Please specify one of \"MGI\"",
                       " or \"MAF\"")
        stop(memo)
    }
    
    # Check if x input is a data frame
    if(!is.data.frame(x))
    {
        memo <- paste0("argument supplied to x is not an object of class",
                       " data frame, attempting to coerce")
        warning(memo)
        x <- as.data.frame(x)
    }
    
    # check for duplicate elements in x
    if(nrow(unique(x)) != nrow(x))
    {
        warning("Detected duplicate rows in x, was this expected?")
    }
    
    # Check if y input is a data frame
    if(!is.null(y))
    {
        # Check y input if data frame
        if(is.data.frame(y))
        {
            if(!all(colnames(y) %in% c('Prop', 'trans_tranv')))
            {
                memo <- paste0("Did not detect correct column names in",
                               "input to y, missing one of \"Prop\",",
                               "\"trans_tranv\"")
                stop(memo)
            }
            
            y$Prop <- as.numeric(as.character(y$Prop))
        }
        
        # Check y input if vector
        if(is.vector(y))
        {
            y <- as.data.frame(y)
            y$trans_tranv <- rownames(y)
            colnames(y) <- c('Prop', 'trans_tranv')
            
            if(typeof(y$Prop) != "double" & typeof(y$Prop) != "numeric")
            {
                memo <- paste0("Values found in y are not of type double",
                               " or numeric!")
                stop(memo)
            }
        }
        
        if(!is.data.frame(y))
        {
            memo <- paste0("input to y is not an object of class data frame",
                           " or named vector")
            stop(memo)
        }
    }
    
    # Check columns of x input and change to internal format
    if(file_type == 'MGI')
    {
        # Check that columns are named appropriatley, if not print error
        proper_names <- c("reference", "variant", "sample")
        if(all(proper_names %in% colnames(x)))
        {
            message("Found appropriate columns")
        } else {
            memo <- paste0("Could not find all columns requested, missing ",
                           "one of \"reference\", \"variant\", \"sample\"")
            stop(memo)
        }
        
        x <- x[,c('reference', 'variant', 'sample')]
        
    } else if(file_type == 'MAF') {
        proper_names <- c("Tumor_Sample_Barcode", "Reference_Allele",
                          "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
        if(all(proper_names %in% colnames(x)))
        {
            message("Found appropriate columns")
        } else {
            memo <- paste0("Could not find all columns requested, missing one ",
                           "of \"Tumor_Sample_Barcode\", \"Reference_Allele\",",
                           " \"Tumor_Seq_Allele1\", \"Tumor_Seq_Allele2\"")
            stop(memo)
        }
        # Convert MAF file to internal format
        x <- TvTi_convMaf(x)
    } else {
        memo <- paste0("TvTi requires a fileType specification, please",
                       "specify one of \"MAF\" or \"MGI\" based on the",
                       "argument supplied to parameter x. See docs for help.")
        stop(memo)
    }
    
    # Remove any indels present in the data
    x <- TvTi_rmIndel(x)
    # Warn about multi nucleotide codes
    x <- TvTi_rmMnuc(x)
    
    # Check that reference and variant columns only contain the proper codes
    ref_codes <- c('A', 'C', 'G', 'T', '-', 0)
    if(!all(toupper(x$reference) %in% toupper(ref_codes)))
    {
        memo <- paste0("Unrecognized Base Detected in reference column, ",
                       "expected values are: ", toString(ref_codes))
        stop(memo)
    } else if(!all(toupper(x$variant) %in% toupper(ref_codes))) {
        memo <- paste0("Unrecognized Base Detected in reference column, ",
                       "expected values are: ", toString(ref_codes))
        stop(memo)
    }
    
    # check y input for proper row names
    if(!is.null(y))
    {
        trans_tranv_names <- c("A->C or T->G (TV)", "A->G or T->C (TI)",
                               "A->T or T->A (TV)", "G->A or C->T (TI)",
                               "G->C or C->G (TV)", "G->T or C->A (TV)")
        if(!all(rownames(y) %in% trans_tranv_names))
        {
            memo <- paste0("Did not detect a value for all combinations of ",
                           "transitions/transversions, please specify input ",
                           "for each of the following levels: ",
                           toString(trans_tranv_names))
            stop(memo)
        }
        
        # check that y sums to 1 (i.e. its a true proportion among all elements)
        if(round(sum(y$Prop), digits=2) != 1)
        {
            stop("The sum of elements in y should equal 1")
        }
    }
    
    # Check input data to clinDat
    if(!is.null(z))
    {
        if(!is.data.frame(z))
        {
            stop("Did not detect a data frame for input to clinDat")
        }
        z <- droplevels(z)
        
        if(!all(c('sample', 'variable', 'value') %in% colnames(z)))
        {
            stop("Did not detect correct sample names in clinDat")
        }
        
        if(!all(levels(x$sample) %in% levels(z$sample)))
        {
            memo <- paste0("Found a sample supplied to clinData not found",
                           " in the data frame supplied to x")
            warning(memo)
        }
    }
    
    return(list('input1'=x, 'input2'=y, 'input3'=z))
}

#' Remove indels
#'
#' Given a data frame with columns reference and variants remove all indels
#' from data
#' @name TvTi_rmIndel
#' @param x Object of class data frame containing columns 'reference', 'variant'
#' @noRd
#' @return Object of class data frame with indels removed

TvTi_rmIndel <- function(x)
{
    original_size <- nrow(x)
    
    # Find and remove insertions and deletions
    x <- x[grep('-|0', x$reference, perl=TRUE, invert=TRUE),]
    x <- x[grep('-|0', x$variant, perl=TRUE, invert=TRUE),]
    
    new_size <- nrow(x)
    
    # Print message if indels have been removed
    if(new_size != original_size)
    {
        message("Removed ", original_size - new_size, " indels present in data")
    }
    
    return(x)
}

#' Remove multinucleotide codes
#'
#' Given a data frame with columns reference and variants remove all
#' multinucleotides from data
#' @name TvTi_rmMnuc
#' @param x Object of class data frame containing columns 'reference', 'variant'
#' @noRd
#' @return Object of class data frame with multi nucleotide codes removed

TvTi_rmMnuc <- function(x)
{
    original_size <- nrow(x)
    
    # Find and remove multi nucleotide codes
    x <- x[grep('[ACGT]{2,}', x$reference, perl=TRUE, invert=TRUE, ignore.case=TRUE),]
    x <- x[grep('[ACGT]{2,}', x$variant, perl=TRUE, invert=TRUE, ignore.case=TRUE),]
    
    new_size <- nrow(x)
    
    if(new_size != original_size)
    {
        memo <- paste0("Multi Nucleotide codes are not currently supported, ", 
                       "removed: ", original_size - new_size,
                       " multi-nucleotides present in the data")
        warning(memo) 
    }
    
    return(x)
}

#' Construct a lolliplot
#'
#' Given a data frame construct a plot displaying mutations on a transcript
#' framework.
#' @name lolliplot
#' @param x Object of class data frame with rows representing mutations. The
#' data frame must contain columns with the following names "transcript_name",
#' "gene", and "amino_acid_change". Values in the "transcript_name" column must
#' represent an ensembl transcript id and values in the "amino_acid_change"
#' column must be in p.notation (see details).
#' @param y Object of class data frame with rows representing mutations. The
#' data frame must contain columns with the following names "transcript_name"
#' and "amino_acid_change". Values in the "transcript_name" column must
#' represent an ensembl transcript id and values in the "amino_acid_change"
#' column must be in p. notation (optional, see details).
#' @param z Object of class data frame with rows representing regions of
#' interest. The data frame must contain columns with the following names
#' "description", "start", "stop" (optional see details).
#' @param fillCol Character string specifying the column name of the argument
#' supplied to parameter x on which to colour the lollis representing mutations
#' (see details).
#' @param labelCol Character string specifying the column name of the argument
#' supplied to parameter x from which to extract and display text corresponding
#' to mutations (see details).
#' @param txtAngle Integer specifying the angle of label text to be plotted if
#' an argument is supplied to the labelCol parameter.
#' @param txtSize Integer specifying the size of label text to be plotted if an
#' argument is supplied to the labelCol parameter.
#' @param pntSize Integer specifying the size of lolli points representing
#' mutations.
#' @param proteinColour Character string specifying the background colour of the
#' protein.
#' @param obsA.rep.fact Numeric value representing the repulsive factor for the
#' lollis plotted, which were derived from the argument supplied to parameter x
#' (see details and vignette).
#' @param obsA.rep.dist.lmt Numberic value representing the repulsive distance
#' limit for the lollis plotted, which were derived from the argument supplied
#' to parameter x (see details and vignette).
#' @param obsA.attr.fact Numeric value representing the attraction factor for
#' the lollis plotted, which were derived from the argument supplied to
#' parameter x (see details and vignette).
#' @param obsA.adj.max Numeric value representing the max position adjustment
#' for the lollis plotted, which were derived from the argument supplied to
#' parameter x (see details and vignette).
#' @param obsA.adj.lmt Numeric value representing the adjustment limit for the
#' lollis plotted, which were derived from the argument supplied to parameter x
#' (see details and vignette).
#' @param obsA.iter.max Integer representing the number of iterations of
#' position adjustments for the lollis plotted, which were derived from the 
#' argument supplied to parameter x (see details and vignette).
#' @param obsB.rep.fact Numeric value representing the repulsive factor for the
#' lollis plotted, which were derived from the argument supplied to parameter y
#' (see details and vignette).
#' @param obsB.rep.dist.lmt Numberic value representing the repulsive distance
#' limit for the lollis plotted, which were derived from the argument supplied
#' to parameter y (see details and vignette).
#' @param obsB.attr.fact Numeric value representing the attraction factor for
#' the lollis plotted, which were derived from the argument supplied to
#' parameter y (see details and vignette).
#' @param obsB.adj.max Numeric value representing the max position adjustment
#' for the lollis plotted, which were derived from the argument supplied to
#' parameter y (see details and vignette).
#' @param obsB.adj.lmt Numeric value representing the adjustment limit for the
#' lollis plotted, which were derived from the argument supplied to parameter y
#' (see details and vignette).
#' @param obsB.iter.max Integer representing the number of iterations of
#' position adjustments for the lollis plotted, which were derived from the 
#' argument supplied to parameter y (see details and vignette).
#' @param sideChain Boolean specifying if amino acid sidechain data should be
#' plotted in lieu of protein domains (see details).
#' @param species A valid species from which to retrieve protein domain and
#' sequence data for a given transcript (see details).
#' @param maxLolliStack Integer specifying the cutoff for the maximum number of
#' lollis allowed to be stacked at a single position.
#' @param plotLayer Valid ggplot2 layer to be added to the plot.
#' @param paletteA Character vector specifying colours for protein domains,
#' valid only if sideChain==FALSE.
#' @param paletteB Character vector specifying colours for lollis representing
#' mutations, valid only if argument is supplied to fillCol.
#' @param host Host to connect to for biomaRt queries (see details).
#' @param out Character vector specifying the the object to output, one of
#' "data", "grob", or "plot", defaults to "plot" (see returns).
#' @details lolliplot is a function designed to display mutation information in
#' the context of a protien identified by an ensembl transcript id. The
#' lolliplot function will query ensembl via biomart to retrieve sequence and
#' domain information in order to construct a representation of a protein and
#' therefore requires an internet connection. A value must be supplied to the
#' species parameter (defaults to hsapiens) in order for a successful biomart
#' query. Valid arguments to this field are those species with datasts available
#' via ensembl. please specify species in lowercase without a period
#' (i.e. hsapiens instead of H.sapiens), lolliplot will inform the user of
#' available species if input to the species parameter is not recognized.
#' Further lolliplot will build a protein framework based on sequence data
#' obtained from biomaRt, by default this will default to the latest ensembl
#' version. In order for the most accurate representation the annotation version
#' of the mutations given to lolliplot should match the annotation version used 
#' by biomaRt. The annotation version used by biomaRt can be changed via the 
#' host paramter (see vignette for more details).
#' 
#' lolliplot is capable of plotting two seperate sets of data on the protein
#' representation specified by parameters `x` and `y`, the data supplied to
#' these parameters will be plotted on the top and bottom of the protein
#' respectively. Note that input to these parameters is expected to correspond
#' to a single ensembl transcript and that values in the "amino_acid_change"
#' columns are required to be in p. notation (i.e. p.V600E). Further lolliplot
#' is able to plot custom domain annotation if supplied via the parameter `z`,
#' this will override domain information obtained from biomart.
#' 
#' lolliplot uses a forcefield model from the package FField to attract and 
#' repulse lollis. The parameters for this force field model are set to
#' reasonable defaults however may be adjusted via the obsA... and obsB...
#' family of parameters. Please see the package FField available on cran for
#' a description of these parameters. Note that the time to construct the
#' lolliplot will in large part depend on the number of mutations and the values
#' supplied to the forcefield parameters.
#' 
#' @examples
#' # Create input data
#' data <- brcaMAF[brcaMAF$Hugo_Symbol == 'TP53',c('Hugo_Symbol', 'amino_acid_change_WU')]
#' data <- as.data.frame(cbind(data, 'ENST00000269305'))
#' colnames(data) <- c('gene', 'amino_acid_change', 'transcript_name')
#'
#' # Call lolliplot
#' lolliplot(data)
#' @return One of the following, a list of dataframes containing data to be
#' plotted, a grob object, or a plot.
#' @export

lolliplot <- function(x, y=NULL, z=NULL, fillCol=NULL, labelCol=NULL,
                      txtAngle=45, txtSize=5, pntSize=4,
                      proteinColour='#999999', obsA.rep.fact=5000,
                      obsA.rep.dist.lmt=500, obsA.attr.fact=.1, obsA.adj.max=.1,
                      obsA.adj.lmt=.5, obsA.iter.max=50000, obsB.rep.fact=5000,
                      obsB.rep.dist.lmt=500, obsB.attr.fact=.1, obsB.adj.max=.1,
                      obsB.adj.lmt=.5, obsB.iter.max=50000,
                      sideChain=FALSE, species="hsapiens",
                      maxLolliStack=NULL, plotLayer=NULL, paletteA=NULL,
                      paletteB=NULL, host="www.ensembl.org", out="plot")
{
    memo <- paste("lolliplot() has been deprecated and may be removed in a future release",
                  "Please use Lolliplot() with a capital L instead!")
    message(memo)
    
    # Perform quality check
    input <- lolliplot_qual(x, y, z)
    x <- input[[1]]
    y <- input[[2]]
    z <- input[[3]]
    
    # extract transcript id and subset data y on that id if it exists
    transcriptID <- as.character(x$transcript_name[1])
    if(!is.null(y))
    {
        y <- y[y$transcript_name == transcriptID,]
    }
    
    # extract HUGO gene name
    gene <- as.character(x$gene[1])
    
    # Obtain length of protein
    result <- lolliplot_transcriptID2codingSeq(transcriptID,
                                               species=species,
                                               host=host)
    codingSeq <- result$coding
    cdsLen <- result$cds_length
    
    
    # Get the sequence length in AA, perform quality checks along the way
    residueSeq <- lolliplot_DNAconv(codingSeq, to="residue")    
    # If it is requested grab the sidechain information and bind to residues
    if(sideChain==TRUE)
    {
        sidechain <- lolliplot_DNAconv(codingSeq, to="sidechain")
        AAsequence <- cbind(sidechain, residueSeq)
        AAsequence <- as.data.frame(AAsequence)
        AAsequence$coord <- seq(from=1, to=nrow(AAsequence))
    } else {
        AAsequence <- NULL
    }
    
    # if there are any stop codons remove them as they are not considered part
    # of the protein
    if(any(residueSeq %in% c("OPAL", "OCHRE", "AMBER")))
    {
        stopRes <- c("OPAL", "OCHRE", "AMBER")
        residueSeq <- residueSeq[-which(residueSeq %in% stopRes)]
        if(!is.null(AAsequence))
        {
            AAsequence <- AAsequence[-which(AAsequence$residueSeq %in% stopRes),]
        }
    }
    
    # grab the length of the protein in Amino Acids
    proteinLength <- length(residueSeq)    
    
    # if z is specified plot that instead of fetching the domain information
    if(!is.null(z))
    {
        geneData <- lolliplot_constructGene(gene, z, proteinLength)
    } else {
        # extract protien domain data
        protein_domain <- lolliplot_fetchDomain(transcriptID,
                                                species=species,
                                                host=host)
        
        # construct gene from data collected
        geneData <- lolliplot_constructGene(gene, protein_domain, proteinLength)
    }
    
    # construct data frame of observed mutations for top track
    observed_mutation <- lolliplot_mutationObs(x, 'top', fillCol, labelCol,
                                               obsA.rep.fact, obsA.rep.dist.lmt,
                                               obsA.attr.fact, obsA.adj.max,
                                               obsA.adj.lmt, obsA.iter.max)
    observed_mutation <- lolliplot_reduceLolli(observed_mutation,
                                               maxLolliStack)
    
    # construct data frame of observed mutations for bottom track
    if(!is.null(y))
    {
        observed_mutation2 <- lolliplot_mutationObs(y, 'bottom', fillCol,
                                                    labelCol, obsB.rep.fact,
                                                    obsB.rep.dist.lmt,
                                                    obsB.attr.fact,
                                                    obsB.adj.max, obsB.adj.lmt,
                                                    obsB.iter.max)
        observed_mutation2 <- lolliplot_reduceLolli(observed_mutation2,
                                                    maxLolliStack)
    } else {
        observed_mutation2 <- NULL
    }
    
    # construct the lolliplot
    plot <- lolliplot_buildMain(geneData, length, observed_mutation,
                                observed_mutation2,fillCol, labelCol,
                                txtAngle, txtSize, pntSize,
                                proteinColour, AAsequence,
                                plot_sidechain=sideChain, layers=plotLayer,
                                paletteA=paletteA, paletteB=paletteB)
    
    # Decide what to output
    dataOut <- list("gene"=geneData,
                    "observed_mutation"=observed_mutation,
                    "observed_mutation2"=observed_mutation2)
    output <- multi_selectOut(data=dataOut, plot=plot, out=out)
    return(output)
}

#' Convert AA to side chain classification
#' 
#' Given the 1 letter code an amino acid, return the side chian classification
#' @name lolliplot_AA2sidechain
#' @param x Character of length 1 giving the 1 letter amino acid code
#' @return Object of class character

lolliplot_AA2sidechain <- function(x)
{
    # Coerce all AA changes to uppercase and then apply switch statement
    x <- toupper(x)
    x <- switch(EXPR=x, "F"="Nonpolar", "L"="Nonpolar", "S"="Polar",
                "Y"="Polar", "C"="Polar", "W"="Nonpolar", "L"="Nonpolar",
                "P"="Nonpolar", "H"="Basic", "Q"="Polar", "R"="Basic",
                "I"="Nonpolar", "M"="Nonpolar", "T"="Polar", "N"="Polar",
                "K"="Basic", "S"="Polar", "R"= "Basic", "V"="Nonpolar",
                "A"="Nonpolar", "D"= "Acidic", "E"="Acidic", "G"="Polar")
    
    return(x)
}

#' Construct Lolliplot
#'
#' Construct Lolliplot given gene and mutation data
#' @name lolliplot_buildMain
#' @param gene_data object of class dataframe giving protien domain and gene
#' information
#' @param length integer specifying the length of the protien in amino acids
#' @param mutation_observed object of class data frame specifying mutations
#' observed in input file
#' @param mutation_observed2 optional object of class data frame specifying
#' additional mutations for bottom track
#' @param fill_value character string specifying the column on which to colour
#' mutation points
#' @param label_column character string specifying the column containing the
#' labels to attach to mutation points
#' @param plot_text_angle numeric value specifying the angle of text to be
#' plotted
#' @param plot_text_size numeric value specifying the size of text to be plotted
#' @param point_size numeric value specigying the size of mutation points
#' @param gene_colour color to shade plotted gene
#' @param sequence_data object of class dataframe giving AA sequence, sidechain,
#' and coord required if plot_sidechain is true
#' @param plot_sidechain boolean specifying whether to plot the AA sidechain
#' instead of domain information
#' @param layers additional ggplot2 layers to plot
#' @param paletteA Character vector specifying colours for gene features
#' @param paletteB Character vector specifying colours for lolli features
#' @return a ggplot2 object
#' @import ggplot2

lolliplot_buildMain <- function(gene_data, length, mutation_observed,
                                mutation_observed2, fill_value, label_column,
                                plot_text_angle, plot_text_size, point_size,
                                gene_colour, sequence_data,
                                plot_sidechain=FALSE,layers=NULL,
                                paletteA=NULL, paletteB=NULL)
{
    # build the various features of the lolliplot
    
    # Build gene base either using domain information
    # or AA sidechain information
    if(plot_sidechain == TRUE)
    {
        sequence_data$coord_start <-
            as.numeric(as.character(sequence_data$coord)) - 1
        sequence_data$coord_end <- as.numeric(as.character(sequence_data$coord))
        gene_plot <- geom_rect(data=sequence_data,
                               mapping=aes_string(xmin='coord_start',
                                                  xmax='coord_end', ymin=-.1,
                                                  ymax=.1, fill='sidechain'))
        domain_plot <- geom_blank()
    } else {
        gene_plot <- geom_rect(data=gene_data[1,],
                               mapping=aes_string(xmin='pos_from',
                                                  xmax='pos_to',
                                                  ymin='height_min',
                                                  ymax='height_max'),
                               fill='#999999', colour='#000000')
        # Take into account there might not be any domains
        if(nrow(gene_data) == 1)
        {
            domain_plot <- geom_blank()
        } else {
            domain_plot <- geom_rect(data=gene_data[-1,],
                                     mapping=aes_string(xmin='pos_from',
                                                        xmax='pos_to',
                                                        ymin='height_min',
                                                        ymax='height_max',
                                                        fill='Domain'),
                                     alpha=1, colour='black')
        }
    }
    
    # Build the Observed track
    observed_plot <- geom_point(data=mutation_observed,
                                mapping=aes_string(x='coord_x_dodge',
                                                   y='coord_y_dodge',
                                                   colour=fill_value),
                                size=point_size)
    observed_line <- geom_segment(data=mutation_observed,
                                  mapping=aes_string(x='mutation_coord',
                                                     y=.1, xend='coord_x_dodge',
                                                     yend=.3))
    observed_line_2 <- geom_segment(data=mutation_observed,
                                    mapping=aes_string(x='coord_x_dodge', y=.3,
                                                       xend='coord_x_dodge',
                                                       yend='coord_y_dodge'))
    
    # Miscelaneous features
    title <- ggtitle(gene_data[1,1])
    x_label <- xlab('Amino Acid Position')
    
    # add a theme and guide to the plot
    theme <- theme(legend.position='bottom',
                   legend.direction='vertical',
                   legend.box='horizontal',
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.title.y=element_blank())
    guide <- guides(colour=guide_legend(ncol=2), fill=guide_legend(ncol=2))
    
    # set colours manually if these are specified
    if(!is.null(paletteA))
    {
        gene_features_fill <- scale_fill_manual(values=paletteA)
    } else {
        gene_features_fill <- geom_blank()
    }
    
    if(!is.null(paletteB))
    {
        lolli_features_fill <- scale_colour_manual(values=paletteB)
    } else {
        lolli_features_fill <- geom_blank()
    }
    
    # construct the plot with or without 2nd observed track
    if(is.null(mutation_observed2))
    {
        y_limits <- ylim(c(-.1, max(mutation_observed$coord_y_dodge) + .1))
        tmp <- data.frame(xmin=0, xmax=0, ymin=0, ymax=0, x=0, y=0, xend=0, yend=0)
        p1 <- ggplot(data=tmp, aes(xmin=0, xmax=0, ymin=0, ymax=0, x=0, y=0, xend=0, yend=0)) +
            gene_plot + domain_plot + observed_line_2 + observed_line +
            observed_plot + x_label + title + y_limits + theme_bw() + theme +
            guide + layers + gene_features_fill + lolli_features_fill
    } else {
        y_limits <- ylim(c(min(mutation_observed2$coord_y_dodge) - .1,
                           max(mutation_observed$coord_y_dodge) + .1))
        if(any(colnames(mutation_observed2) %in% fill_value))
        {
            observed2_plot <- geom_point(data=mutation_observed2,
                                         mapping=aes_string(x='coord_x_dodge',
                                                            y='coord_y_dodge',
                                                            colour=fill_value),
                                         size=point_size)
        } else {
            observed2_plot <- geom_point(data=mutation_observed2,
                                         mapping=aes_string(x='coord_x_dodge',
                                                            y='coord_y_dodge'),
                                         size=point_size)
        }
        observed2_line <- geom_segment(data=mutation_observed2,
                                       mapping=aes_string(x='mutation_coord',
                                                          y=-.1,
                                                          xend='coord_x_dodge',
                                                          yend=-.3))
        observed2_line_2 <- geom_segment(data=mutation_observed2,
                                         mapping=aes_string(x='coord_x_dodge',
                                                            y=-.3,
                                                            xend='coord_x_dodge',
                                                            yend='coord_y_dodge'))
        
        tmp <- data.frame(xmin=0, xmax=0, ymin=0, ymax=0, x=0, y=0, xend=0, yend=0)
        p1 <- ggplot(tmp, aes(xmin=0, xmax=0, ymin=0, ymax=0, x=0, y=0, xend=0, yend=0)) +
            gene_plot + domain_plot + observed_line + observed_line_2 +
            observed_plot + observed2_line + observed2_line_2 + observed2_plot +
            x_label + title + y_limits + theme_bw() + theme + guide + layers +
            gene_features_fill + lolli_features_fill
    }
    
    # If a label column is specified plot labels
    if(any(colnames(mutation_observed) %in% "labels"))
    {
        mutation_observed$y_label_offset <-
            mutation_observed$coord_y_dodge + .01
        p1 <- p1 + geom_text(data=mutation_observed,
                             mapping=aes_string(x='coord_x_dodge',
                                                y='y_label_offset',
                                                label='labels'),
                             angle=plot_text_angle, size=plot_text_size,
                             vjust=1, hjust=0)
    }
    if(any(colnames(mutation_observed2) %in% "labels"))
    {
        mutation_observed2$y_label_offset <-
            mutation_observed2$coord_y_dodge - .01
        p1 <- p1 + geom_text(data=mutation_observed2,
                             mapping=aes_string(x='coord_x_dodge',
                                                y='y_label_offset',
                                                label='labels'),
                             angle=plot_text_angle, size=plot_text_size,
                             vjust=0, hjust=1)
    }
    
    return(p1)
}

#' Convert Codon to AA
#' 
#' Convert a Codon to the appropriate amino acid
#' @name lolliplot_Codon2AA
#' @param x Character string of length 1 giving the DNA codon to convert
#' @return Character corresponding to the residue for the given codon

lolliplot_Codon2AA <- function(x)
{
    # Convert codons to single AA code
    x <- toupper(x)
    x <- switch(x, "TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L", "CTT"="L",
                "CTC"="L", "CTA"="L", "CTG"="L", "ATT"="I", "ATC"="I",
                "ATA"="I", "ATG"="M", "GTT"="V", "GTC"="V", "GTA"="V",
                "GTG"="V", "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S",
                "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P", "ACT"="T",
                "ACC"="T", "ACA"="T", "ACG"="T", "GCT"="A", "GCC"="A",
                "GCA"="A", "GCG"="A", "TAT"="Y", "TAC"="Y", "TAA"="OCHRE",
                "TAG"="AMBER", "CAT"="H", "CAC"="H", "CAA"="Q", "CAG"="Q",
                "AAT"="N", "AAC"="N", "AAA"="K", "AAG"="K", "GAT"="D",
                "GAC"="D", "GAA"="E", "GAG"="E", "TGT"="C", "TGC"="C",
                "TGA"="OPAL", "TGG"="W", "CGT"="R", "CGC"="R", "CGA"="R",
                "CGG"="R", "AGT"="S", "AGC"="S", "AGA"="R", "AGG"="R",
                "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G")
    
    return(x)
}

#' Construct gene information
#' 
#' Build gene for input into lolliplot_buildMain
#' @name lolliplot_constructGene
#' @param gene character string specifying gene name
#' @param domain_data object of class data frame specifying protien domain
#' information, obtained from lolliplot_fetchDomain, should contain columns
#' giving "description", "start", "end"
#' @param length integer specifying length of transcript in amino acids
#' @return object of class data frame giving gene and domain information

lolliplot_constructGene <- function(gene, domain_data, length)
{
    # message
    message("Constructing gene track")
    
    # Construct basic gene information, if there are no domains return the gene
    gene <- data.frame(Domain=gene, pos_from=1, pos_to=length, nest=1)
    if(nrow(na.omit(domain_data)) == 0)
    {
        gene$height_min <- .1/(as.numeric(gene$nest))
        gene$height_max <- -.1/(as.numeric(gene$nest))
        gene$pos_from <- as.numeric(gene$pos_from)
        gene$pos_to <- as.numeric(gene$pos_to)
        return(gene)
    }
    
    # rename columns for domain_data and make sure description column is
    # not a factor
    colnames(domain_data) <- c("description", "start", "end")
    domain_data$description <- as.character(domain_data$description)
    
    # quality check of domain data
    if(max(domain_data$end) > length)
    {
        memo <- paste0("The end position of a domain: ",  max(domain_data$end),
                       " is exceeding the length of the protein:", length)
        warning(memo)
    } else if(min(domain_data$start) < 1) {
        memo <- paste0("The start position of a domain:",
                       min(domain_data$start),
                       "is less than the start of the protein", 1)
        warning(memo)
    }
    
    # Check that start coordinates are always less than the end coordinates
    if(any(domain_data$start >= domain_data$end))
    {
        memo <- paste0("Found a start position greater than an end position",
                       " in the protein features track. Check input to Z or",
                       "results of the biomaRt query using dataOut==TRUE.")
        warning(memo)
    }
    
    # determine which regions are overlapping and annotate which nest domain is
    # sort on start
    domain_data$start <- as.numeric(domain_data$start)
    domain_data$end <- as.numeric(domain_data$end)
    domain_data <- domain_data[order(domain_data$start),]
    
    # annotate nests
    nest <- vector('numeric')
    end <- vector('numeric')
    for(i in 1:nrow(domain_data))
    {
        # Remove from end any values <= gene$start[i]
        idx <- domain_data$start[i] < end
        end <- end[idx]
        
        nest <- c(nest, length(end))
        end <- c(end, domain_data$end[i])
    }
    
    # add this nest information to the data frame
    domain_data$nest <- nest + 1
    colnames(domain_data) <- c("Domain", "pos_from", "pos_to", "nest")
    
    # combine gene and domain information
    gene <- rbind(gene, domain_data)
    
    # annotate display heights based on nesting and make sure coord are numeric
    gene$height_min <- .1/(as.numeric(gene$nest))
    gene$height_max <- -.1/(as.numeric(gene$nest))
    gene$pos_from <- as.numeric(gene$pos_from)
    gene$pos_to <- as.numeric(gene$pos_to)
    
    return(gene)
}

#' Convert DNA character string
#' 
#' Convert a character string of nucleotides to amino acids or side chain class
#' @name lolliplot_DNAconv
#' @param x Character string of nucleotides to convert
#' @param to Character string specifying conversion to do, one of "codon", 
#' "residue", "sidechain"
#' @return Converted string of nucleotides as character vector

lolliplot_DNAconv <- function(x, to="residue")
{
    # check if given character string is a multiple of 3
    if(nchar(x)%%3 != 0)
    {
        memo <- paste0("Coding sequence retrieved for given ensembl transctipt",
                       ", is not a multiple of three. output may not be,",
                       " accurate!")
        warning(memo)
    }
    
    # split the character string into codons
    codon <- substring(x, seq(1,nchar(x), 3), seq(3, nchar(x), 3))
    if(toupper(to)=="CODON")
    {
        return(as.character(codon))
    }
    
    # convert the codons into amino acid residues
    residue <- sapply(codon, lolliplot_Codon2AA)
    if(toupper(to)=="RESIDUE")
    {
        return(as.character(residue))
    }
    
    # convert the residues into sidechain classifications
    sidechain <- sapply(residue, lolliplot_AA2sidechain)
    if(toupper(to)=="SIDECHAIN")
    {
        return(as.character(sidechain))
    }
    
    # return a warning if code gets this far
    memo <- paste0("did not recognize input to variable \"to\",",
                   " returning residue data")
    warning(memo)
    return(as.character(residue))
}

#' dodge coordinates
#'
#' given amino acid position dodge on x axis
#' @name lolliplot_dodgeCoordX
#' @param x numeric vector of position coordinates on x axis
#' @param rep.fact repulsive factor for plotted mutations observed track
#' @param rep.dist.lmt repulsive distance limit for plotted mutations observed
#' track
#' @param attr.fact attraction factor for plotted mutations observed track
#' @param adj.max maximum position change for each iteration observed track
#' @param adj.lmt position adjustment limit which simulation stops observed
#' track
#' @param iter.max maximum iterations beyond which to stop the simulation
#' observed track
#' @return numeric vector of dodged position coordinates on x axis
#' @importFrom FField FFieldPtRep

lolliplot_dodgeCoordX <- function(x, rep.fact=5000, rep.dist.lmt=500,
                                  attr.fact=.1, adj.max=.1, adj.lmt=.5,
                                  iter.max=50000)
{
    # Format into data frame with columns as x and y
    x <- as.data.frame(cbind(x, 0))
    colnames(x) <- c('x', 'y')
    
    # Forcefield does not work with only 1 point, check for this first
    if(nrow(x) < 2)
    {
        return(x$x)
    }
    
    # take the data frame and apply a repulsive force to coordinates
    x <- FField::FFieldPtRep(x, rep.fact=rep.fact, rep.dist.lmt=rep.dist.lmt,
                             attr.fact=attr.fact, adj.max=adj.max,
                             adj.lmt=adj.lmt, iter.max=iter.max)
    
    return(x$x)
}

#' dodge coordinates
#'
#' given a data frame, dodge x coordinates ontop of each other
#' @name lolliplot_dodgeCoordY
#' @param x data frame containing columns coord_x_dodge
#' @param track character vector, one of "top", "bottom" specifying whether to
#' dodge in a positive or negative fashion
#' @return numeric vector of dodged position coordinates on y axis

lolliplot_dodgeCoordY <- function(x, track='top')
{
    for(i in 1:length(x$coord_x_dodge))
    {
        if(track == 'top' & i == 1)
        {
            pos <- .3
            orig_pos <- .3
            pos_change <- .1
        } else if(track == 'bottom' & i == 1) {
            pos <- -.3
            orig_pos <- -.3
            pos_change <- -.1
        }
        
        if(i == 1)
        {
            y_axis_vec <- c(pos)
            next
        } else {
            x_coord_a <- x$coord_x_dodge[i-1]
            x_coord_b <- x$coord_x_dodge[i]
        }
        
        if(x_coord_b == x_coord_a)
        {
            new_pos <- pos + pos_change
            y_axis_pos <- new_pos
            y_axis_vec <- c(y_axis_vec, y_axis_pos)
            pos <- new_pos
            next
        } else {
            pos <- orig_pos
            y_axis_vec <- c(y_axis_vec, pos)
        }
    }
    
    return(y_axis_vec)
}

#' fetch protein domains
#' 
#' Retrieve protein domains given ensembl transcript ID
#' @name lolliplot_fetchDomain
#' @param transcriptID String specifying ensembl transcript id
#' @param species character string to use when searching for ensemblMart dataset
#' @param host Host to connect to.
#' @return data frame of protien domains and start/stop coordinates
#' @importFrom biomaRt useMart
#' @importFrom biomaRt listDatasets
#' @importFrom biomaRt useDataset
#' @importFrom biomaRt getBM
#' 
lolliplot_fetchDomain <- function(transcriptID,
                                  species="hsapiens",
                                  host="www.ensembl.org")
{
    # display message
    message("Querying biomaRt for protein domains")
    
    # Load in mart
    ensembl_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                                     host="www.ensembl.org")
    
    # select proper data set given regexp print warnings if unexpected out occur
    dataset <- biomaRt::listDatasets(ensembl_mart)$dataset
    index <- which(grepl(species, dataset))
    if(length(index)>1)
    {
        memo <- paste0(species, " Matches more than one dataset for the",
                       " ensembl mart, please specify a species in the, ",
                       "following format: hsapiens")
        stop(memo)
    } else if(length(index)==0) {
        memo <- paste0(species, " does not appear to be supported by biomaRt",
                       "if you beleive this to be in error please modify", 
                       "you're input to to conform to this format: hsapiens")
        stop(memo)
    }
    ensembl_mart <- biomaRt::useDataset(as.character(dataset[index]),
                                        mart=ensembl_mart)
    
    # Apply various filters using vector of values
    filters <- c("ensembl_transcript_id")
    values <- as.list(c(transcriptID))
    
    # Select attributes to retrieve (protein domain, start, stop)
    attributes <- c("interpro_description",
                    "interpro_start",
                    "interpro_end")
    
    # Retrieve data
    result <- biomaRt::getBM(attributes=attributes, filters=filters,
                             values=values, mart=ensembl_mart)
    
    return(result)
}

#' format mutation observations
#'
#' Create a data frame of mutation observations
#' @name lolliplot_mutationObs
#' @param x object of class data frame with columns trv_type and amino acid
#' change
#' @param track character string specifying one to 'top', 'bottom' to specify
#' proper track
#' @param fill_value character string giving the name of the column to shade
#' variants on
#' @param label_column character string specifying column containing text
#' information to be plotted
#' @param rep.fact repulsive factor for plotted mutations observed track
#' @param rep.dist.lmt repulsive distance limit for plotted mutations observed
#' track
#' @param attr.fact attraction factor for plotted mutations observed track
#' @param adj.max maximum position change for each iteration observed track
#' @param adj.lmt position adjustment limit which simulation stops observed
#' track
#' @param iter.max maximum iterations beyond which to stop the simulation
#' observed track
#' @return object of class data frame giving mutation observations

lolliplot_mutationObs <- function(x, track, fill_value, label_column,
                                  rep.fact, rep.dist.lmt, attr.fact, adj.max,
                                  adj.lmt, iter.max)
{
    # Remove variants within an intronic or splice site region
    if(any(grepl("^e", x$amino_acid_change, ignore.case=TRUE, perl=TRUE)))
    {   
        # save original data frame size before subset for message
        origDim <- nrow(x)
        
        # remove regions with AA change starting with e (i.e. intronic/splice)
        x <- x[-which(grepl("^e", x$amino_acid_change,
                            ignore.case=TRUE, perl=TRUE)),]
        
        newDim <- nrow(x)
        
        # print update message
        memo <- paste0("Removed ", origDim - newDim,
                       " variants not within a residue")
        message(memo)
        
        # if the removal has removed all rows print an error
        if(newDim == 0)
        {
            memo <- paste0("Did not detect any residues, please check input", 
                           " lolliplot must have at least one valid residue",
                           " present!")
            stop(memo)
        }
    }
    
    # extract the mutation types and set a flag specifying they are present
    if(any(colnames(x) %in% fill_value))
    {
        fill_value_flag <- TRUE
        fill <- as.character(x[,eval(fill_value)])
    } else {
        fill_value_flag <- FALSE
    }
    
    # extract the mutation coordinates
    mutation_coord <- x$amino_acid_change
    if(all(grepl("p\\.", mutation_coord)))
    {
        message("Detected p. notation for amino_acid_change")
        mutation_coord <- as.numeric(gsub("p\\.[*a-zA-z]*(\\d+).*?$", "\\1",
                                          mutation_coord, perl=TRUE))
    } else if(all(grepl("c\\.", mutation_coord))) {
        memo <- paste0("c. notation is not currently supported",
                       " please specify amino acid change in p. notation")
        stop(memo)
    } else {
        memo <- paste0("Could not determine notation type for ",
                       "column \"amino_acid_change\", please check input.", 
                       "Expecting p. notation: ex. p.R383A")
        stop(memo)
    }
    
    # combine mutation type and mutation coord into a data frame
    if(fill_value_flag)
    {
        mutation_data <- as.data.frame(cbind(mutation_coord, fill))
        colnames(mutation_data) <- c('mutation_coord', eval(fill_value))
    } else {
        mutation_data <- as.data.frame(mutation_coord)
        colnames(mutation_data) <- c('mutation_coord')
    }
    mutation_data$mutation_coord <-
        as.numeric(as.character(mutation_data$mutation_coord))
    
    # add extra column giving height of Y axis for points to be plotted
    if(track == 'top')
    {
        mutation_data$height_max <- .3
    } else if (track == 'bottom') {
        mutation_data$height_min <- -.3
    } else {
        stop("Fatal error: incorrect track type specified")
    }
    
    # extract the mutation types and set a flag specifying they are present
    if(any(colnames(x) %in% label_column))
    {
        label_column_flag <- TRUE
        mutation_data$labels <- as.character(x[,eval(label_column)])
    } else {
        label_column_flag <- FALSE
    }
    
    # Dodge mutation coordinates on the x axis
    if(track == 'top')
    {
        memo <- paste0("applying force field to observed mutations for",
                       " top track. This will take time if n is large",
                       ", see vignette for tips")
        message(memo)
    } else if (track == 'bottom') {
        memo <- paste0("applying force field to observed mutations for",
                       " bottom track. This will take time if n is large",
                       ", see vignette for tips")        
        message(memo)
    }
    mutation_data <- mutation_data[order(mutation_coord),] 
    mutation_data$coord_x_dodge <- 
        lolliplot_dodgeCoordX(as.vector(mutation_data$mutation_coord),
                              rep.fact=rep.fact,
                              rep.dist.lmt=rep.dist.lmt,
                              attr.fact=attr.fact,
                              adj.max=adj.max,
                              adj.lmt=adj.lmt,
                              iter.max=iter.max)
    
    # Dodge y coordinates
    if(track == 'top')
    {
        mutation_data$coord_y_dodge <- lolliplot_dodgeCoordY(mutation_data,
                                                             track='top')
    } else if(track == 'bottom') {
        mutation_data$coord_y_dodge <- lolliplot_dodgeCoordY(mutation_data,
                                                             track='bottom')
    }
    
    return(mutation_data)
}

#' Check input to lolliplot
#'
#' Perform Basic quality checks for lolliplot input
#' @name lolliplot_qual
#' @param x object of class data frame containing columns transcript_name, gene,
#' and amino_acid_change and rows denoting mutations
#' @param y object of class data frame containing columns transcript_name, and
#' amino_acid_change and rows denoting mutations
#' @param z Object of class data frame containing columns "description", "start",
#'  "stop" specifying gene regions to highlight
#' @return objects passing basic quality checks

lolliplot_qual <- function(x, y, z)
{
    # Check input to x
    if(!is.data.frame(x))
    {
        message("Input to x is not a data frame, attempting to coerce")
        x <- as.data.frame(x)
        x <- droplevels(x)
    }
    
    # Check for correct columns in x
    if(!all(c('transcript_name', 'gene', 'amino_acid_change') %in% colnames(x)))
    {
        stop("Did not detect correct columns in x,
             missing one of transcript_name, gene, amino_acid_change")
    }
    
    # Make sure columns in x are of the proper class
    x$transcript_name <- as.factor(x$transcript_name)
    x$gene <- as.factor(x$gene)
    x$amino_acid_change <- as.factor(x$amino_acid_change)
    
    # Check that "transcript_name" in x contains only 1 transcript
    if(length(unique(x$transcript_name)) != 1)
    {
        stop("Detected more than 1 transcript in input to x")
    }
    
    # Check input to y
    if(!is.null(y))
    {
        # is y a data frame?
        if(!is.data.frame(y))
        {
            message(y, "is not a data frame, attempting to coerce")
            y <- as.data.frame(y)
            y <- droplevels(y)
        }
        
        # does y have correct columns?
        if(!all(c('transcript_name', 'amino_acid_change') %in% colnames(y)))
        {
            stop("Did not detect correct columns in y, missing one of
                 transcript_name, amino_acid_change")
        }
        
        # make sure columns in y are of proper class
        y$transcript_name <- as.factor(y$transcript_name)
        y$amino_acid_change <- as.factor(y$amino_acid_change)
        }
    
    # Check input to z
    if(!is.null(z))
    {
        # is z a data frame?
        if(!is.data.frame(z))
        {
            memo <- paste0("Input to z is not a data frame",
                           ", attempting to coerce")
            warning(memo)
            z <- as.data.frame(z)
            z <- droplevels(z)
        }
        
        # does z contain correct columns
        if(!all(c("description", "start", "stop") %in% colnames(z)))
        {
            memo <- paste0("Did not detect correct columns in input to z, ",
                           "missing one of \"description\", \"start\",",
                           " \"stop\"")
            stop(memo)
        }
        
        # make sure column class of z are of proper type
        z$description <- as.factor(z$description)
        z$start <- as.numeric(as.character(z$start))
        z$stop <- as.numeric(as.character(z$stop))
    }
    
    return(list(x, y, z))
    }

#' Reduce Lolli
#' 
#' Reduce lollis stacked ontop of each other to the amount specified
#' @name lolliplot_reduceLolli
#' @param x Data frame with column name mutation_coord to reduce lollis on
#' @param max Integer specifying the maximum number of lollis to allow
#' @return Object of class data frame taking the reduced form of x
#' @importFrom plyr count

lolliplot_reduceLolli <- function(x, max=NULL)
{
    # if max is null no reduction is required
    if(is.null(max))
    {
        return(x)
    }
    
    # Get a frequency of counts from x
    coordFreq <- plyr::count(x$mutation_coord)
    colnames(coordFreq) <- c("x", "freq")
    
    keep <- vector('numeric')
    # Loop through mutations keeping only what is under max
    for(i in 1:nrow(coordFreq))
    {
        index <- which(x$mutation_coord == coordFreq$x[i])
        index <- index[1:max]
        keep <- c(keep, index)
    }
    
    # subset the input on what we want to keep
    x <- x[keep,]
    
    return(x)
}

#' fetch protein length
#' 
#' Retrieve protein length from ensembl database given enseml transcript id
#' @name lolliplot_transcriptID2codingSeq
#' @param transcriptID character string giving ensembl transcript id
#' @param species character string to use when searching for ensemblMart dataset
#' @param host Host to connect to.
#' @return length in residues of ensembl transcript id
#' @importFrom biomaRt useMart
#' @importFrom biomaRt listDatasets
#' @importFrom biomaRt useDataset
#' @importFrom biomaRt getBM

lolliplot_transcriptID2codingSeq <- function(transcriptID,
                                             species="hsapiens",
                                             host="www.ensembl.org")
{
    # display mesage
    memo <- paste0("Using the following host: ", host, " for biomaRt queries",
                   " to change the ensembl annotation version alter this",
                   " parameter!")
    message(memo)
    message("Querying biomaRt for transcript sequence")
    
    # Load in mart
    ensembl_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                                     host=host)
    
    # select proper data set given regexp print warnings if unexpected out occur
    dataset <- biomaRt::listDatasets(ensembl_mart)$dataset
    index <- which(grepl(species, dataset))
    if(length(index)>1)
    {
        memo <- paste0(species, " Matches more than one dataset for the",
                       " ensembl mart, please specify a species in the, ",
                       "following format: hsapiens")
        stop(memo)
    } else if(length(index)==0) {
        valid_species <- toString(gsub("_gene_ensembl",
                                       "",
                                       dataset))
        
        memo <- paste0(species, " does not appear to be supported by biomaRt",
                       " please specify one of the following species:",
                       valid_species)
        stop(memo)
    }
    ensembl_mart <- biomaRt::useDataset(as.character(dataset[index]),
                                        mart=ensembl_mart)
    
    # Apply various filters using vector of values
    filters <- c("ensembl_transcript_id")
    ensg_id <- as.character(transcriptID)
    
    # Select attributes to retrieve coding dna sequence
    attributes <- as.list(c("coding","cds_length"))
    
    # Retrieve data
    result <- biomaRt::getBM(attributes=attributes, filters=filters,
                             values=ensg_id, mart=ensembl_mart)
    
    return(as.list(result))
}

