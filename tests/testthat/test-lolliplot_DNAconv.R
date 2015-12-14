test_that("lolliplot_DNAconv identifies if an incomplete codon is given in the sequence", {
    x <- c("AAATT")
    expect_warning(lolliplot_DNAconv(x), "not a multiple of three")
})

test_that("lolliplot_DNAconv splits the input character string into codon lengths", {
    x <- c("AAATTT")
    expect_equal(length(lolliplot_DNAconv(x)), 2)
})

test_that("lolliplot_DNAconv converts condons into residues", {
    x <- c("AAA")
    expect_equal(lolliplot_DNAconv(x), "K")
})

test_that("lolliplot_DNAconv converts residues into sidechains", {
    x <- c("AAA")
    expect_equal(lolliplot_DNAconv(x, to="sidechain"), "Basic")
})

test_that("lolliplot_DNAconv throws an warning if input to the to parameter is unrecognized", {
    x <- c("AAA")
    expect_warning(lolliplot_DNAconv(x, to="incorrect"), "did not recognize")
})