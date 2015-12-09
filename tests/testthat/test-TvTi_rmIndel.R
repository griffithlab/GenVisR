test_that("TvTi_rmIndel correctly identifies and removes indels", {
    x <- data.frame("reference"=c('-', 'A'), "variant"=c('A', 'A'))    
    expect_equal(nrow(suppressWarnings(TvTi_rmIndel(x))), 1)
    
    x <- data.frame("reference"=c('A', 'A'), "variant"=c('0', 'A'))    
    expect_equal(nrow(suppressWarnings(TvTi_rmIndel(x))), 1)    
})

test_that("TvTi_rmIndel correctly warns if indels are detected", {
    x <- data.frame("reference"=c('-', 'A'), "variant"=c('A', 'A'))
    
    expect_message(TvTi_rmIndel(x), "indels present")
})