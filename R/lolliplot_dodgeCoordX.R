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
