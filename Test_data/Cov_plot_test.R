require("GenomicRanges")
setwd("/Users/awagner/Workspace/R/GGgenome/Test_data")

# need a biostrings object for reference
require(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19

# need transcript data for reference
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# need Granges object 
gr <- GRanges(seqnames=c("chr16"), ranges=IRanges(start=c(67063050), end=c(67134958)), strand=strand(c("+")))
#gr <- GRanges(seqnames=c("chr16"), ranges=IRanges(start=c(67020842), end=c(67136898)), strand=strand(c("+")))

# need coverage data and extra stuff to make plot work
cov <- read.delim('CBFB_TCGA_Cov.bed')
RNA <- apply(cov[,4:ncol(cov)], 1, mean)
cov$cov <- RNA
cov <- as.data.frame(cbind(cov[,3], cov[,'cov']))
colnames(cov) <- c('end', 'cov')
#cov1 <- cov[,c(3,4)]
#colnames(cov1) <- c("end", "cov")
#cov2 <- cov[,c(3,5)]
#colnames(cov2) <- c("end", "cov")
data <- list("RNA" = cov)
cov2 <- cov[7480:7500,]
data2 <- list("RNA" = cov2)


#test <- gene_plot(txdb, gr, genome, reduce=F, transformIntronic=T)
#plot_coverage(data, txdb, gr, genome, reduce=T)
#plot_track('CBFB'=gene, 'Cov'=test_plot, 'Cov'=test_plot, border='green', axis_align='width', width_ratio=c(1,10))

# test <- function(x)
# {
#   x$start <- x$end
#   return(x)
# }
# cov2 <- test(cov2)
# 
# master <- gene_plot(txdb, gr, genome, reduce=FALSE, transformIntronic=TRUE, output_transInt_table=TRUE)
# test <- adply(cov2, 1, map_coord_space, master=master, .parallel=FALSE)

test <- plot_coverage(data2, txdb, gr, genome, reduce=F, transformIntronic=T)
