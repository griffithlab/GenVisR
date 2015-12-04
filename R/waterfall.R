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
#' @param mutBurden an optional data frame containing columns sample, mut_burden
#' with sample levels matching those supplied in x
#' @param main.recurrence_cutoff integer specifying removal of entries not seen
#' in at least "x" percent of samples
#' @param main.grid a boolean value to overlay a grid on the primary plot
#' @param main.label_x a boolean value to plot samples on the x axis
#' @param main.gene_label_size an integer specifying the size of labels on Y
#' axis
#' @param coverageSpace an integer specifying the size in bp of the genome
#' covered from which mutations could be called
#' @param file_type a character string specifying the file format of the data
#' frame, one of "MGI", "MAF", "Custom"
#' @param main.genes a character vector specifying genes to plot
#' @param drop_mutation Boolean specifying whether to drop unused
#' "mutation type" levels from the legend
#' @param rmv_silent Boolean specifying wheter to remove silent mutations from
#' the left side and main plot
#' @param main.label_col Character string specifying an optional column name
#' from which to derive cell labels from
#' @param main.plot_label_size Integer specifying the size of labels for cells
#' in the main panel
#' @param main.plot_label_angle Integer specifying the degree of rotation for
#' text if main.label_col is specified
#' @param main.palette Character vector specifying colors to fill on mutation
#' type
#' @param sampRecur.layers Additional ggplot2 layers to plot on the sample
#' recurence chart
#' @param main.layers Additional ggplot2 layers to plot on the main panel
#' @param mutRecur.layers Additional ggplot2 layers to plot on the mutation
#' burden data
#' @param variant_class_order character vector giving the hierarchical order of
#' mutation types to plot (required if file_type="Custom")
#' @param main.samples character vector containing sample names to subset the
#' input from x on
#' @examples
#' waterfall(brcaMAF)
#' @return a grob for plotting
#' @export

waterfall <- function(x, clinData=NULL, clinLegCol=1, clinVarCol=NULL,
                      clinVarOrder=NULL, clinLayer=NULL, mutBurden=NULL,
                      main.recurrence_cutoff=0, main.grid=TRUE,
                      main.label_x=FALSE, main.gene_label_size=8,
                      coverageSpace=44100000, file_type='MAF', main.genes=NULL,
                      main.samples=NULL, drop_mutation=FALSE, rmv_silent=FALSE,
                      main.label_col=NULL, main.plot_label_size=4,
                      main.palette=NULL, sampRecur.layers=NULL,
                      main.layers=NULL, mutRecur.layers=NULL,
                      main.plot_label_angle=0, variant_class_order=NULL)
{
    # Perform data quality checks and conversions
    inputDat <- waterfall_qual(x, clinData, mutBurden, file_type=file_type,
                               label_col=main.label_col)
    data_frame <- inputDat[[1]]
    clinData <- inputDat[[2]]
    mutBurden <- inputDat[[3]]

    # Set flag if it is desirable to plot cell text
    if(!is.null(main.label_col))
    {
        main.plot_label_flag <- TRUE
    } else {
        main.plot_label_flag <- FALSE
    }

    # If it is requested subset the input data on a sample list
    if(!is.null(main.samples))
    {
        data_frame <- waterfall_sampAlt(data_frame, main.samples)
    }

    # add in a count of mutations at the sample level before anything is
    # stripped out and save for mutation recurrence plot
    data_frame2 <- waterfall_calcMutFreq(data_frame[,c('sample', 'gene',
                                                       'trv_type')])

    # Subset the data to remove silent mutations if specified
    if(rmv_silent==TRUE)
    {
        data_frame <- waterfall_rmvSilent(data_frame)
    }

    # Subset the data based on a vector of genes if supplied
    if(!is.null(main.genes))
    {
        data_frame <- waterfall_geneAlt(data_frame, main.genes)
    }

    # Remove trv_type that are not the most deleterious for a given gene/sample
    data_frame <- waterfall_hierarchyTRV(data_frame, file_type,
                                         variant_class_order)

    # Subset the data based on the recurrence of mutations at the gene level
    data_frame <- waterfall_geneRecurCutoff(data_frame, main.recurrence_cutoff)

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
                 either: the samples given in x, or main.samples")
        }

        mutBurden$sample <- factor(mutBurden$sample, levels=sample_order)
        p3 <- waterfall_buildMutBurden_B(mutBurden, layers=mutRecur.layers)
    } else {
        data_frame2$sample <- factor(data_frame2$sample,
                                     levels=sample_order)
        p3 <- waterfall_buildMutBurden_A(data_frame2, coverageSpace,
                                         layers=mutRecur.layers)
    }

    # Plot the Left Bar Chart
    p2 <- waterfall_buildGenePrevelance(data_frame, layers=sampRecur.layers)

    # if there are any NA values in the data frame for a gene, give these NA
    # values a gene name so they are plotted properly
    data_frame <- waterfall_NA2gene(data_frame)

    # Plot the Heatmap
    if(is.null(clinData))
    {
        p1 <- waterfall_buildMain(data_frame, grid=main.grid,
                                  label_x=main.label_x,
                                  gene_label_size=main.gene_label_size,
                                  file_type=file_type,
                                  drop_mutation=drop_mutation,
                                  plot_x_title=TRUE,
                                  plot_label=main.plot_label_flag,
                                  plot_label_size=main.plot_label_size,
                                  plot_palette=main.palette, layers=main.layers,
                                  plot_label_angle=main.plot_label_angle)
    } else if(!is.null(clinData)) {
        p1 <- waterfall_buildMain(data_frame, grid=main.grid,
                                  label_x=main.label_x,
                                  gene_label_size=main.gene_label_size,
                                  file_type=file_type,
                                  drop_mutation=drop_mutation,
                                  plot_x_title=FALSE,
                                  plot_label=main.plot_label_flag,
                                  plot_label_size=main.plot_label_size,
                                  plot_palette=main.palette, layers=main.layers,
                                  plot_label_angle=main.plot_label_angle)
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
