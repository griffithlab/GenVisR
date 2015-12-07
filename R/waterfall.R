#' Construct a waterfall plot
#'
#' Given a data frame construct a water fall plot showing the mutation burden
#' and mutation type on a gene and sample level.
#' @name waterfall
#' @param x Object of class data frame representing annotated mutations. The 
#' data frame supplied must have one of the following sets of column names
#' ("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant Classification") for
#' fileType="MAF", ("sample","gene_name","trv_type") for fileType="MGI" or
#' ("sample", "gene", "variant_class") for fileType="Custom". This columns 
#' should represent samples in a cohort, gene with mutation, and the mutation
#' type respectively.
#' @param clinData Object of class data frame with rows representing clinical
#' data. The data frame should be in "long format" and columns must be names as
#' "sample", "variable", and "value" (optional see details and vignette).
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
#' @param mutBurden Object of class data frame containing columns "sample",
#' "mut_burden" with sample levels matching those supplied in x.
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
#' @param coverageSpace Integer specifying the size in bp of the genome
#' covered by sequence data from which mutations could be called
#' (see details and vignette).
#' @param fileType Character string specifying the file format of the data
#' frame specified to parameter `x`, one of "MGI", "MAF", "Custom" 
#' (see details and vignette).
#' @param plotGenes Character vector specifying genes to plot. If not null genes
#' not specified within this character vector are removed.
#' @param mainDropMut Boolean specifying whether to drop unused
#' "mutation type" levels from the legend.
#' @param rmvSilent Boolean specifying if silent mutations should be removed 
#' from the plot.
#' @param mainLabelCol Character string specifying a column name from the 
#' argument supplied to parameter `x` from which to derive cell labels from
#' (see details and vignette).
#' @param mainLabelSize Integer specifying the size of text labels for cells
#' in the main plot. Valid only if argument is supplied to the parameter
#' `mainLabelCol`.
#' @param mainLabelAngle Integer specifying the degree of rotation for
#' text labels. Valid only if argument is supplied to the parameter
#' `mainLabelCol`.
#' @param mainPalette Character vector specifying colours for mutation types
#' plotted in the main plot, must specify a colour for each mutation type
#' plotted.
#' @param sampRecurLayer Valid ggplot2 layer to be added to the left sub-plot.
#' @param mainLayer Valid ggplot2 layer to be added to the main plot.
#' @param mutBurdenLayer Valid ggplot2 layer to be added to the top sub-plot.
#' @param variant_class_order Character vector specifying the hierarchical order
#' of mutation types to plot, required if file_type == "Custom"
#' (see details and vignette).
#' @param plotSamples Character vector specifying samples to plot. If not null
#' all other samples not specified within this parameter are removed.
#' @param dataOut Boolean specifying whether to output the data to be passed to
#' ggplot instead of plotting.
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
#' waterfall(brcaMAF)
#' @return A graphic object.
#' @export

waterfall <- function(x, clinData=NULL, clinLegCol=1, clinVarCol=NULL,
                      clinVarOrder=NULL, clinLayer=NULL, mutBurden=NULL,
                      mainRecurCutoff=0, mainGrid=TRUE,
                      mainXlabel=FALSE, main_geneLabSize=8,
                      coverageSpace=44100000, fileType='MAF', plotGenes=NULL,
                      plotSamples=NULL, mainDropMut=FALSE, rmvSilent=FALSE,
                      mainLabelCol=NULL, mainLabelSize=4,
                      mainPalette=NULL, sampRecurLayer=NULL,
                      mainLayer=NULL, mutBurdenLayer=NULL,
                      mainLabelAngle=0, variant_class_order=NULL, dataOut=FALSE)
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

    # reorder the genes based on frequency of mutations in the gene
    gene_sorted <- waterfall_geneSort(data_frame)
    data_frame$gene <- factor(data_frame$gene, levels=gene_sorted)

    # reorder the samples based on hiearchial sort on ordered gene list
    sample_order <- waterfall_sampSort(data_frame)
    data_frame$sample <- factor(data_frame$sample, levels=sample_order)

    # Reorder the sample levels in data_frame2 to match the main plot's levels,
    # and then plot the top margin plot
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

    # Plot the Left Bar Chart
    p2 <- waterfall_buildGenePrevelance(data_frame, layers=sampRecurLayer)

    # if there are any NA values in the data frame for a gene, give these NA
    # values a gene name so they are plotted properly
    data_frame <- waterfall_NA2gene(data_frame)
    
    if(isTrue(dataOut))
    {
        return(data_frame)
    }

    # Plot the Heatmap
    if(is.null(clinData))
    {
        p1 <- waterfall_buildMain(data_frame, grid=mainGrid,
                                  label_x=mainXlabel,
                                  gene_label_size=main_geneLabSize,
                                  file_type=fileType,
                                  drop_mutation=mainDropMut,
                                  plot_x_title=TRUE,
                                  plot_label=main.plot_label_flag,
                                  plot_label_size=mainLabelSize,
                                  plot_palette=mainPalette, layers=mainLayer,
                                  plot_label_angle=mainLabelAngle)
    } else if(!is.null(clinData)) {
        p1 <- waterfall_buildMain(data_frame, grid=mainGrid,
                                  label_x=mainXlabel,
                                  gene_label_size=main_geneLabSize,
                                  file_type=fileType,
                                  drop_mutation=mainDropMut,
                                  plot_x_title=FALSE,
                                  plot_label=main.plot_label_flag,
                                  plot_label_size=mainLabelSize,
                                  plot_palette=mainPalette, layers=mainLayer,
                                  plot_label_angle=mainLabelAngle)
    }

    # Plot any clinical data if it is specified
    if(!is.null(clinData))
    {
        # match the levels of sample in y to conform to the main plot
        clinData$sample <- factor(clinData$sample, levels=sample_order)
        p4 <- multi_buildClin(clinData, clin.legend.col=clinLegCol, 
                              clin.var.colour=clinVarCol, 
                              clin.var.order=clinVarOrder, 
                              clin.layers=clinLayer)

        # Align all plots and return as 1 plot
        pA <- waterfall_align(p2, p1, p3, p4)
        return(grid::grid.draw(pA))
    }

    # Align the Plots and return as 1 plot
    pA <- waterfall_align(p2, p1, p3)

    return(grid::grid.draw(pA))
}
