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
#' @param proportions ggplot object displaying proportion of mutation types "optional"
#' @param section_heights Heights of each section (should sum to one)
#' @return a grob object
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid nullGrob

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
    
    if (!is.null(proportions)) {
        prop_grob <- ggplot2::ggplotGrob(proportions)
        ind_legend <- grep("guide", prop_grob$layout$name)
        prop_legend <- prop_grob[["grobs"]][[ind_legend]]
        prop_width <- sum(prop_legend$width)
        prop_grob <- ggplot2::ggplotGrob(proportions + theme(legend.position="none"))
        
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
        if (!is.null(proportions)) {
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
        if (!is.null(proportions)) {
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
        if (!is.null(proportions)) {
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
        if (!is.null(proportions)) {
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
#' 
#' @description Builds a ggplot object showing individuals' mutational profile
#' 
#' @return a ggplot object
#' 
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