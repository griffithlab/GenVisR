#' Create Region Table
#' 
#' Create a master region table by merging isoforms
#' @name mergeRegions
#' @param txdb A TxDb object for a genome
#' @param gr Granges object specifying the ROI
#' @return Master region table data frame

mergeRegions <- function(gene_features, gr, base=exp(1)){
  #extract preserved intron ranges, exons, and intron buffer regions, and update Type
  master <- as.data.frame(do.call("rbind", gene_features))[,c('start','end','width','Type')]
  
  #combine Intron (compressed) spaces
  master <- unique(master[order(master$start, master$end),])
  master$Type <- as.character(master$Type)

  gr.start <- as.data.frame(ranges(gr))[1,1]
  gr.stop <- as.data.frame(ranges(gr))[1,2]

  i <- 1
  while (i < nrow(master)){
    r.start <- master[i,2] + 1
    r.stop <- master[i+1,1] - 1
    width <- r.stop - r.start + 1
    if (width < 1){
      new.start <- master[i,1]
      new.stop <- master[i+1,2]
      new.width <- new.stop - new.start + 1
      master[i,1:3] <- c(new.start,new.stop,new.width)
      master[i,4] <- 'Merged'
      master <- master[-(i+1),]
    }else{
      i <- i + 1
    }
  }
  
  if (gr.start < master[1,'start']){
    m.start <- master[1,'start']
    row <- master[1,]
    row[1:3] <- c(gr.start, m.start - 1, log(m.start-gr.start, base=base))
    row[4] <- 'Intron'
    master <- rbind(master, row)
  }
  for (i in 1:(nrow(master)-1)){
    r.start <- master[i,2] + 1
    r.stop <- master[i+1, 1] - 1
    width <- r.stop - r.start + 1
    row <- master[1,]
    row[1:3] <- c(r.start, r.stop, log(width, base=base))
    row[4] <- 'Intron'
    master <- rbind(master, row)
  }
  if (gr.stop > master[i+1, 2]){
    r.start <- master[i+1, 2] + 1
    row <- master[1,]
    row[1:3] <- c(r.start, gr.stop, log(gr.stop - r.start, base=base))
    row[4] <- 'Intron'
    master <- rbind(master, row)
  }
  
  master <- master[order(master$start, master$end),]
  
  # Transform original coordinates into new space and bind to master data frame
  for(i in 1:nrow(master))
  {
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
  
  return(master)
}