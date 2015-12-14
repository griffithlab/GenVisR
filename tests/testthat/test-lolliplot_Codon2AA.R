test_that("lolliplot_Codon2AA Converts codons to amino acids", {
    x <- c("TTG")
    out <- lolliplot_Codon2AA(x)
    expect_equal(out, "L")
    
    x <- c("ttg")
    out <- lolliplot_Codon2AA(x)
    expect_equal(out, "L")
})

test_that("lolliplot_Codon2AA outputs NA if a codon is not recognized", {
    x <- c("ZZZ")
    out <- lolliplot_Codon2AA(x)
    expect_equal(out, NULL)
})