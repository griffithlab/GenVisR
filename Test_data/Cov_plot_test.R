#setwd("/Users/zskidmor/GGgenome/Test_data")
setwd("/Users/awagner/Workspace/R/GGgenome/Test_data")
# need a biostrings object for reference
require(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19

# need transcript data for reference
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# need Granges object 
gr <- GRanges(seqnames=c("chr16"), ranges=IRanges(start=c(67063051), end=c(67134782)), strand=strand(c("+")))
#gr <- GRanges(seqnames=c("chr16"), ranges=IRanges(start=c(67063051), end=c(67065051)), strand=strand(c("+")))

# need coverage data and extra stuff to make plot work
cov <- read.delim('CBFB_TCGA_Cov.bed')
RNA <- apply(cov[,4:ncol(cov)], 1, mean)
cov$cov <- RNA
cov <- as.data.frame(cbind(cov[,3], cov[,'cov']))
colnames(cov) <- c('end', 'cov')
data <- list("RNA" = cov)

# test <- gene_plot(txdb, gr, genome, reduce=T, transform=c('Intron','UTR','CDS'))
test <- plot_coverage(data, txdb, gr, genome, reduce=F, transform=c('Intron','UTR','CDS'))
# test <- plot_track('CBFB'=gene, 'Cov'=test_plot, 'Cov'=test_plot, border='green', axis_align='width', width_ratio=c(1,10))
