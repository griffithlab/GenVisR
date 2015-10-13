#' Merge Typed Region Tables
#'
#' Create a master region table by merging isoforms
#' @name mergeTypes
#' @param master An unsorted dataframe of CDS and UTR elements before merging
#' @return master A sorted dataframe of merged CDS and UTR elements

mergeTypes <- function(master){
    master <- master[order(master$start, master$end),]
    if (nrow(master) >= 2){# If there are fewer than 2 rows, then there aren't any conflicts to resolve.
        for (i in 1:(nrow(master)-1)){
            first.end <- master[i,2]
            second.start <- master[i+1,1]
            dist <- second.start - first.end - 1
            if (dist < 0){
                second.isInternal <- master[i+1,2] < master[i,2]
                if (master[i,]$Type == 'CDS'){
                    if(second.isInternal){# Here `second` is a UTR completely encompassed by a CDS (`first`)
                        master<-master[-(i+1),]# Remove the internal UTR region from the table
                        return(mergeTypes(master))
                    }
                    master[i+1,1] <- first.end + 1
                    master[i+1,]$width <- master[i+1, 2] - master[i+1, 1] + 1
                }else{
                    master[i, 2] <- second.start - 1
                    master[i,]$width <- master[i, 2] - master[i, 1] + 1
                    if(second.isInternal){# Here `second` is a CDS completely encompassed by a UTR (`first`)
                        UTR.endSegment <- master[i,]
                        UTR.endSegment$start <- master[i+1, 2] + 1
                        UTR.endSegment$end <- first.end# Unchanged from beginning of loop
                        UTR.endSegment$width <- UTR.endSegment$end - UTR.endSegment$start + 1
                        master <- rbind(master, UTR.endSegment)# Add the right UTR segment to dataframe
                        return(mergeTypes(master))
                    }
                }
            }
        }
    }
    return(master)
}
