test_that("multi_cytobandRet retrieves data in the expected format", {
    out <- suppressWarnings(multi_cytobandRet("hg19"))
    expect <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
    expect_equal(colnames(out), expect)
})