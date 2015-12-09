test_that("TvTi_rmMnuc correctly identifies and removes multi-nucleotide codes", {
    x <- data.frame("reference"=c("AA", "T"), "variant"=c("AT", "C"))
    expect_equal(nrow(suppressWarnings(TvTi_rmMnuc(x))), 1)

    x <- data.frame("reference"=c("gg", "t"), "variant"=c("cc", "g"))
    expect_equal(nrow(suppressWarnings(TvTi_rmMnuc(x))), 1)
})

test_that("TvTi_rmMnuc warns if it finds multi-nucleotide codes", {
    x <- data.frame("reference"=c("AA", "T"), "variant"=c("AT", "C"))
    
    expect_warning(TvTi_rmMnuc(x), "multi-nucleotides present")   
})