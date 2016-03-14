test_that("waterfall_geneSort properly orders genes", {
    x <- data.frame(sample=rep('a', 10), gene=c(rep(c('b', 'c', 'd'), 3), 'e'), trv_type=rep('missense', 10))
    out <- waterfall_geneSort(x) 
    expec <- c('e', 'b', 'c', 'd')
    expect_equal(out, expec)
})

test_that("waterfall_geneSort defers to input to geneOrder if specified by the user", {
    x <- data.frame(sample=rep('a', 10), gene=c(rep(c('b', 'c', 'd'), 3), 'e'), trv_type=rep('missense', 10))
    gene_order <- c("d", "b", "c", "e")
    out <- waterfall_geneSort(x, geneOrder=gene_order)
    expec <- rev(gene_order)
    expect_equal(out, expec)
    
    gene_order <- c("none", "b", "d")
    expec <- c("d", "b")
    out <- suppressWarnings(waterfall_geneSort(x, geneOrder=gene_order))
    expect_equal(out, expec)
})

test_that("waterfall_geneSort warns if it finds genes supplied to geneOrder not in the primary input", {
    # test case where no genes in geneOrder are in primary input
    x <- data.frame(sample=rep('a', 10), gene=c(rep(c('b', 'c', 'd'), 3), 'e'), trv_type=rep('missense', 10))
    gene_order <- c("none")
    expect_warning(waterfall_geneSort(x, geneOrder=gene_order), "Did not find any genes")
    
    # test case where some but not all genes in geneOrder are in primary input
    x <- data.frame(sample=rep('a', 10), gene=c(rep(c('b', 'c', 'd'), 3), 'e'), trv_type=rep('missense', 10))
    gene_order <- c("none", "b", "d")
    expect_warning(waterfall_geneSort(x, geneOrder=gene_order), "The following genes were not found")
})