#' Construct a oncoprint
#'
#' This function has been removed, please use Waterfall() (capital W) Tutorial can be found at: https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpz1.252
#' @name waterfall
#' @export

waterfall <- function()
{
    message("This function has been replaced, please use Waterfall() instead, a full tutorial can be found at: https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpz1.252")
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
        
        if(!all(unique(x$sample) %in% unique(z$sample)))
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
#' This function has been removed, please use Lolliplot() (capital L) instead!
#' @name lolliplot
#' @export

lolliplot <- function()
{
    message("This function has been removed due to the removal of the FField package on cran, please use Lolliplot() instead")
}
