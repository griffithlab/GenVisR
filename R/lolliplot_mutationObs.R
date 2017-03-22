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
