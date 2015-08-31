#' Plot a mutation landscape
#' 
#' Plot a mutation landscape plot for a cohort in an annotation file
#' @name mutspec
#' @param x a data frame in annotation format
#' @param clinDat an optional data frame in "long" format giving additional
#' information to be plotted, requires columns "sample", "variable", and "value"
#' @param clin.legend.col an integer specifying the number of columns to plot in
#' the clinical data legend
#' @param clin.var.colour a named character vector specifying the mapping
#' between colors and variables in the clinical data
#' @param clin.var.order a character vector of variables to order the clinical
#' legend by
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
#' frame, one of "MGI", "MAF"
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
#' @param clin.layers Additional ggplot2 layers to plot on the clinical data 
#' @param main.layers Additional ggplot2 layers to plot on the main panel
#' @param mutRecur.layers Additional ggplot2 layers to plot on the mutation
#' burden data
#' @examples
#' mutSpec(brcaMAF)
#' @return a grob for plotting
#' @export

mutSpec <- function(x, clinDat=NULL, clin.legend.col=1, clin.var.colour=NULL,
                    clin.var.order=NULL, mutBurden=NULL,
                    main.recurrence_cutoff=0, main.grid=TRUE,
                    main.label_x=FALSE, main.gene_label_size=8, 
                    coverageSpace=44100000, file_type='MAF', main.genes=NULL,
                    main.samples=NULL, drop_mutation=FALSE, rmv_silent=FALSE,
                    main.label_col=NULL, main.plot_label_size=4,
                    main.palette=NULL, sampRecur.layers=NULL, clin.layers=NULL,
                    main.layers=NULL, mutRecur.layers=NULL,
                    main.plot_label_angle=0)
{ 
    # Perform data quality checks and coversions
    inputDat <- mutSpec.qual(x, clinDat, mutBurden, file_type=file_type,
                             label_col=main.label_col)
    data_frame <- inputDat[[1]]
    clinDat <- inputDat[[2]]
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
        data_frame <- mutSpec.mutation_sample_subset(data_frame, main.samples)
    }
    
    # add in a count of mutations at the sample level before anything is
    # stripped out and save for mutation recurrence plot
    data_frame2 <- mutSpec.add_mutation_counts(data_frame[,c('sample',
                                                             'gene',
                                                             'trv_type')])
    
    # Subset the data to remove silent mutations if specified
    if(rmv_silent==TRUE)
    {
        data_frame <- mutSpec.mutation_silent_rmv(data_frame)
    }
    
    # Remove trv_type that are not the most deleterious for a given gene/sample
    data_frame <- mutSpec.hiearchial_remove_trv_type(data_frame,
                                                     file_type=file_type)
    
    # reorder the genes based on frequency of mutations in the gene
    gene_sorted <- mutSpec.gene_sort(data_frame)
    data_frame$gene <- factor(data_frame$gene, levels=gene_sorted)
    
    # reorder the samples based on hiearchial sort on ordered gene list
    sample_order <- mutSpec.sample_sort(data_frame)
    data_frame$sample <- factor(data_frame$sample, levels=sample_order)
    
    # Subset the data based on the recurrence of mutations at the gene level
    data_frame <- mutSpec.mutation_recurrence_subset(data_frame,
                                                     main.recurrence_cutoff)
    
    # Subset the data based on a vector of genes if supplied
    if(!is.null(main.genes))
    {
        data_frame <- mutSpec.mutation_gene_subset(data_frame, main.genes)
    }
    
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
        p3 <- build.mutRecurB.mutSpec(mutBurden, layers=mutRecur.layers)
        } else {
            data_frame2$sample <- factor(data_frame2$sample,
                                         levels=sample_order)
            p3 <- build.mutRecurA.mutSpec(data_frame2, coverageSpace,
                                          layers=mutRecur.layers)
        }
    
    # Plot the Left Bar Chart
    p2 <- build.mutOccur.mutSpec(data_frame, layers=sampRecur.layers)
    
    # if there are any NA values in the data frame for a gene, give these NA
    # values a gene name so they are plotted properly
    data_frame <- mutSpec.add_gene_to_NA(data_frame)
    
    # Plot the Heatmap
    if(is.null(clinDat))
    {
        p1 <- build.main.mutSpec(data_frame, grid=main.grid,
                                 label_x=main.label_x,
                                 gene_label_size=main.gene_label_size,
                                 file_type=file_type,
                                 drop_mutation=drop_mutation,
                                 plot_x_title=TRUE,
                                 plot_label=main.plot_label_flag,
                                 plot_label_size=main.plot_label_size,
                                 plot_palette=main.palette, layers=main.layers,
                                 plot_label_angle=main.plot_label_angle)
    } else if(!is.null(clinDat)) {
        p1 <- build.main.mutSpec(data_frame, grid=main.grid,
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
    if(!is.null(clinDat))
    {
        # match the levels of sample in y to conform to the main plot
        clinDat$sample <- factor(clinDat$sample, levels=sample_order)
        p4 <- build.clin.mutSpec(clinDat, clin.legend.col=clin.legend.col, 
                                 clin.var.colour=clin.var.colour, 
                                 clin.var.order=clin.var.order,
                                 clin.layers=clin.layers)
        
        # Align all plots and return as 1 plot
        pA <- mutSpec.align_waterfall(p2, p1, p3, p4)
        return(grid::grid.draw(pA))
    }
    
    # Align the Plots and return as 1 plot
    pA <- mutSpec.align_waterfall(p2, p1, p3)
    
    return(grid::grid.draw(pA))
}
