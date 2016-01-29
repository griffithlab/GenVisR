test_that("waterfall_geneAlt checks for proper input", {
    x <- data.frame(sample='a', gene='b', trv_type='missense')
    genes <- c(1) 
    expect_warning(waterfall_geneAlt(x, genes), "not a character vector")
})

test_that("waterfall_geneAlt checks that genes are available for removal", {
    x <- data.frame(sample='a', gene='b', trv_type='missense')
    genes <- c('z')
    expect_warning(waterfall_geneAlt(x, genes), "element not found in x")
})

test_that("waterfall_geneAlt successfully removes genes from input", {
    x <- data.frame(sample=c('a','a'), gene=c('b','c'), trv_type=rep('missense',2))
    genes <- c('b')
    out <- waterfall_geneAlt(x, genes)
    expect_equal(nrow(out), 1)
})
