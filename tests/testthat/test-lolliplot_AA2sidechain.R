test_that("lolliplot_AA2sidechain, converts amino acid codes to sidechain information", {
    x <- c("f")
    out <- lolliplot_AA2sidechain(x)
    expect_equal(out, "Nonpolar")
    
    x <- c("F")
    out <- lolliplot_AA2sidechain(x)
    expect_equal(out, "Nonpolar")
})