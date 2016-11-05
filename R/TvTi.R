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
