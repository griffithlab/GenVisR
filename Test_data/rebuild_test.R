genome <- BSgenome.Hsapiens.UCSC.hg19
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

gr <- GRanges(seqnames=c("chr16"), ranges=IRanges(start=c(67063051), end=c(67134782)), strand=strand(c("+")))
#gr <- GRanges(seqnames=c("chr16"), ranges=IRanges(start=c(67063051), end=c(67065051)), strand=strand(c("+")))
#gr <- GRanges(seqnames=c("chrX"), ranges=IRanges(start=c(73820651), end=c(73852753)), strand=strand(c("-"))) #XIST

test <- plot_coverage(data, txdb, gr, genome, reduce=F, transform=c('Intron','UTR','CDS'))