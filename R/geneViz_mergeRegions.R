#' Create Region Table
#'
#' Create a master region table by merging isoforms
#' @name geneViz_mergeRegions
#' @param gene_features A dataframe specifying features of a gene
#' @param gr Granges object specifying the ROI
#' @param base A vector of log bases to transform the data, corresponding to the
#' elements of transform
#' @param transform A vector of strings designating what objects to log
#' transform
#' @return Master region table data frame
#' @importFrom IRanges ranges

geneViz_mergeRegions <- function(gene_features, gr, base, transform)
{
    #base / transform check
    if(length(base) != length(transform))
    {
        stop("Base vector shorter than transform vector.")
    }

    #extract preserved intron ranges, exons, and intron buffer regions,
    # and update Type
    master <- as.data.frame(do.call("rbind", gene_features))[,c('start','end',
                                                                'width','Type')]

    #order regions and remove duplicates
    master <- unique(master[order(master$start, master$end),])
    master$Type <- as.character(master$Type)

    #merge CDS
    CDS.master <- master[master$Type == 'CDS',]
    CDS.master <- geneViz_mergeTypeRegions(CDS.master)

    #merge UTR
    UTR.master <- master[master$Type == 'UTR',]
    UTR.master <- geneViz_mergeTypeRegions(UTR.master)

    #combine CDS and UTR. CDS regions take precedence if overlap.
    master <- rbind(CDS.master, UTR.master)
    width.init <- master$width
    master <- cbind(master,width.init)
    master <- geneViz_mergeTypes(master)

    #create introns / intergenic space TODO: classify intronic vs. intergenic space
    gr.start <- as.data.frame(IRanges::ranges(gr))[1,1]
    gr.stop <- as.data.frame(IRanges::ranges(gr))[1,2]

    #internal introns / intergenic
    n <- nrow(master)
    if (n >= 2)
    {
        for (i in 1:(n-1))
        {
            r.start <- master[i,2] + 1
            r.stop <- master[i+1, 1] - 1
            width <- r.stop - r.start + 1
            if(width < 1){next}
            row <- master[1,]
            row[c(1:3,5)] <- c(r.start, r.stop, width, width)
            row[4] <- 'Intron'
            master <- rbind(master, row)
        }
    }
    
    #upstream intron / intergenic
    if (gr.start < master[1,'start'])
    {
        m.start <- master[1,'start']
        row <- master[1,]
        width <- m.start-gr.start
        row[c(1:3,5)] <- c(gr.start, m.start - 1, width, width)
        row[4] <- 'Intron'
        master <- rbind(master, row)
    }
    
    #downstream intron / intergenic
    if (gr.stop > master[n, 2])
    {
        r.start <- master[n, 2]
        row <- master[1,]
        width <- gr.stop - r.start
        row[c(1:3,5)] <- c(r.start + 1, gr.stop, width, width)
        row[4] <- 'Intron'
        master <- rbind(master, row)
    }

    master <- master[order(master$start, master$end),]

    # Transform original coordinates into new space and bind
    # to master data frame
    for(i in 1:nrow(master))
    {
        if(master[i,4] %in% transform)
        {
            idx = match(master[i,4], transform)
            master[i,3] <- log(master[i,5], base=base[idx]) * master[i,3]/master[i,5]
        }
        
        if(i == 1)
        {
            trans_start <- 1
            trans_end <- master[1,c('width')] + trans_start
            trans_start_vec <- trans_start
            trans_end_vec <- trans_end
        } else {
            trans_start <- trans_end_vec[i-1]
            trans_end <- master[i,c('width')] + trans_end_vec[i-1]
            trans_start_vec <- c(trans_start_vec, trans_start)
            trans_end_vec <- c(trans_end_vec, trans_end)
        }
    }

    master$trans_start <- trans_start_vec
    master$trans_end <- trans_end_vec
    r <- master$width / master$width.init
    master['c.ratio'] <- max(r)/r
    
    return(master)
}
