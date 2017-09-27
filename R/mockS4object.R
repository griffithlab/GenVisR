#' create mock of S4 object
#' 
#' @name getPosition
#' @rdname getPosition-methods
#' @param object Object of class S4
#' @param ... additional arguments to passed to fill in S4 slots
#' @details function originally from an Rhelp group:
#' http://r.789695.n4.nabble.com/Using-a-mock-of-an-S4-class-td4728547.html and
#' proposed by Martin Morgan.
#' @noRd
mockS4object <- function(class, ..., where=topenv(parent.frame())) { 
    obj <- .Call( 
        methods:::C_new_object, 
        getClassDef(class, where=where) 
    ) 
    
    args = list(...) 
    for (nm in names(args)) 
        slot(obj, nm) = args[[nm]] 
    
    return(obj) 
} 