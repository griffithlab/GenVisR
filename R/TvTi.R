#' plot transitions/transversions
#'
#' Given a data frame with columns reference, variant, and sample construct a
#' transition/transversion plot
#' @name TvTi
#' @param x Object of class data frame containing columns 'sample', reference',
#' 'variant'
#' @param y named vector containing expected transition/transversion proportions
#' with names "A->C or T->G (TV)", "A->G or T->C (TI)", "A->T or T->A (TV)",
#'  "G->A or C->T (TI)", "G->C or C->G (TV)", "G->T or C->A (TV)" or a data
#' frame with column names "Prop", "trans_tranv" and levels of trans_tranv
#'  matching "A->C or T->G (TV)", "A->G or T->C (TI)", "A->T or T->A (TV)",
#'  "G->A or C->T (TI)", "G->C or C->G (TV)", "G->T or C->A (TV)"
#' @param type Object of class character specifying whether to plot the
#' Proportion or Frequency, one of "Proportion", "Frequency"
#' @param label_x_axis boolean specifying wheter to label x axis
#' @param x_axis_text_angle Integer specifying the angle of labels on the x_axis
#' @param palette Character vector of length 6 specifying colors for each
#' trans/tranv type
#' @param file_type Character string specifying the format the input is in,
#' one of 'MAF', 'MGI'
#' @param tvti.layers Additional ggplot2 layers to add to the main panel plot
#' @param expec.layers Additional ggplot2 layers to add to the expected
#' values plot
#' @param sort character specifying one of "sample", "tvti", "none" for
#' plotting samples based on a sample name, transition/transversion freq, or no
#' order
#' @param dataOut Boolean Specifying whether to output the data to be passed
#'  instead of plotting it
#' @param clin.legend.col an integer specifying the number of columns to plot in
#' the clinical data legend
#' @param clin.var.colour a named character vector specifying the mapping
#' between colors and variables in the clinical data
#' @param clin.var.order a character vector of variables to order the clinical
#' legend by
#' @param clin.layer Additional ggplot2 layer to plot on the clinical data
#' @examples
#' TvTi(brcaMAF, type='Frequency',
#' palette=c("#77C55D", "#A461B4", "#C1524B", "#93B5BB", "#4F433F", "#BFA753"),
#' x_axis_text_angle=60)
#' @return ggplot2 object or data frame if dataOut is set to TRUE
#' @export


TvTi <- function(x, y=NULL, clinData=NULL, type='Proportion', label_x_axis=TRUE,
                 x_axis_text_angle=45,
                 palette=c('#D53E4F', '#FC8D59', '#FEE08B', '#E6F598',
                           '#99D594', '#3288BD'),
                 file_type='MAF', tvti.layers=NULL, expec.layers=NULL,
                 sort='none', dataOut=FALSE, clin.legend.col=NULL,
                 clin.var.colour=NULL, clin.var.order=NULL, clin.layer=NULL)
{ 
    # Perform quality checks
    out <- TvTi_qual(x, y, clinData, file_type=file_type)
    x <- out$input1
    y <- out$input2
    clinData <- out$input3

    # add transition/transversion info
    message("annotating transitions and transversions")
    x <- plyr::adply(x, 1, TvTi_annoTransTranv, .progress='text')

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
    } else if(toupper(sort) == toupper('none')){
        sample_order <- levels(x$sample)
    } else {
        memo <- paset0(sort, " is not a valid parameter for sort, please",
                       " specify one of \"sample\", \"tvti\", \"none\"")
        stop(memo)
    }

    # if requested output the data instead of ploting
    if(isTRUE(dataOut))
    {
        return(x)
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
                             x_axis_text_angle=x_axis_text_angle,
                             palette=palette, label_x_axis=label_x_axis,
                             tvti.layers=tvti.layers, expec.layers=NULL,
                             title_x_axis=TRUE)    
    }
   
    
    # Plot Clinical Data if Speccified and build modified TvTi main plot
    if(!is.null(clinData))
    {
        clinData$sample <- factor(clinData$sample, levels=sample_order)
        p3 <- multi_buildClin(clinData, clin.legend.col=clin.legend.col, 
                              clin.var.colour=clin.var.colour, 
                              clin.var.order=clin.var.order,
                              clin.layers=clin.layer)
        
        # Build the Transition/Transversion Plot
        p1 <- TvTi_buildMain(x, y, type=type,
                             x_axis_text_angle=x_axis_text_angle,
                             palette=palette, label_x_axis=label_x_axis,
                             tvti.layers=tvti.layers, expec.layers=NULL,
                             title_x_axis=FALSE)
    } else {
        p3 <- NULL
    }

    if(!is.null(y))
    {
        # If y is input plot the expected values
        p2 <- TvTi_buildMain(y, y, type=type,
                             x_axis_text_angle=x_axis_text_angle,
                             palette=palette, label_x_axis=label_x_axis,
                             plot_expected=TRUE, tvti.layers=NULL,
                             expec.layers=expec.layers)
    } else {
        p2 <- NULL
    }

    finalPlot <- TvTi_alignPlot(p1=p1, p2=p2, p3=p3)
    return(grid::grid.draw(finalPlot))
}
