#' Check input to mutSpec
#' 
#' Perform a data quality check on input to mutSpec
#' @name mutSpec.qual
#' @param x a data frame in annotation format
#' @param y a data frame containing clinical data or a null object
#' @param z a data frame containing mutation burden information or a null object
#' @param file_type Character string specifying the input format to expect in x
#' @param label_col Character string specifying the column name of a label column
#' @return a list of data frames passing quality checks

mutSpec.qual <- function(x, y, z, file_type, label_col)
{
  # Check input data to x
  if(!is.data.frame(x))
  {
    warning("Did not detect a data frame for x, attempting to coerce")
    x <- as.data.frame(x)
  }
  
  # drop unused levels in x
  x <- droplevels(x)
  
  # Convert file type to internal format
  if(toupper(file_type) == toupper("MAF"))
  {
    x <- MAF_to_anno(x, label_col)
  } else if(toupper(file_type) == toupper("MGI"))
  {
    x <- MGI_to_anno(x, label_col)
  } else {
    stop("Unrecognized file_type: ", file_type)
  }
  
  # Check input data to clinDat
  if(!is.null(y))
  {
    y <- droplevels(y)
    
    if(!is.data.frame(y))
    {
      warning("Did not detect a data frame for clinDat, attempting to coerce")
      y <- as.data.frame(y)
    }
    
    if(!all(c('sample', 'variable', 'value') %in% colnames(y)))
    {
      stop("Did not detect correct sample names in clinDat")
    }
  }
  
  # check input data to mutBurden
  if(!is.null(z))
  {
    z <- droplevels(z)
    
    if(!is.data.frame(z))
    {
      warning("Did not detect a data frame for mutBurden, attempting to coerce")
      z <- as.data.frame(z)
    }
    
    if(!all(c('sample', 'mut_burden') %in% colnames(z)))
    {
      stop("Did not detect correct sample names in mutBurden")
    }
  }

  return(list(x, y, z))
}