test_that("waterfall_sampAlt checks for proper class", {
    x <- data.frame(sample='1', gene='egfr', trv_type='missense')
    samples <- c(1)
    expect_warning(waterfall_sampAlt(x, samples), "not a character vector")
})

test_that("waterfall_sampAlt correctly adds in samples not in the original data", {
    x <- data.frame(sample='samp1', gene='egfr', trv_type='missense')
    samples <- c('samp1', 'samp2')
    out <- waterfall_sampAlt(x, samples)
    expec <- data.frame(sample=c('samp2', 'samp1'), gene=c(NA, 'egfr'), trv_type=c(NA, 'missense'))
    expect_equal(out, expec)
})

test_that("waterfall_sampAlt correctly removes samples in the original data", {
    x <- data.frame(sample=c('samp1', 'samp2'), gene=c('egfr', 'tp53'), trv_type=c('missense', 'nonsense'))
    samples <- c('samp1')
    out <- waterfall_sampAlt(x, samples)
    expec <- data.frame(sample='samp1', gene='egfr', trv_type='missense')
    expect_equal(out, expec)   
})