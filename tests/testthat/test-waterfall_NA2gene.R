test_that("waterfall_NA2gene properly assigns NA values in the input to the top gene", {
    x <- data.frame(sample=c('a', 'b', NA), gene=c('c', 'c', NA), trv_type=c('missense', 'missense', NA))
    out <- waterfall_NA2gene(x)
    expect_false(any(is.na(out$gene)))
})