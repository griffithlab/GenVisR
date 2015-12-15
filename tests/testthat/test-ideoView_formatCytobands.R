test_that("ideoView_formatCytobands properly formats input", {
    x <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == "hg19",]
    out <- head(ideoView_formatCytobands(x, "chr1"), 2)
    expect <- data.frame(chrom=c("chr1", "chr1"), chromStart=c(0, 2300000),
                         chromEnd=c(2300000, 5400000), name=c("p36.33", "p36.32"),
                         gieStain=c("gneg", "gpos25"), genome=c("hg19", "hg19"),
                         height_min=c(-0.5, -0.5), height_max=c(0.5, 0.5),
                         alternate=c("top", "bottom"), band_center=c(1150000, 3850000),
                         text_y=c(0.7, -0.7), arm=c("p", "p"),
                         stringsAsFactors=FALSE)
    expect_equivalent(out, expect)
})