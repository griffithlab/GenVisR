test_that("TvTi correctly levels transitions and transversion for alignment with sub-plots", {
    x <- data.frame(sample=c("samp1", "samp1"),
                    reference=c("a", "a"),
                    variant=c("c", "c"))
    out <- suppressWarnings(TvTi(x,
                               y=c("A->C or T->G (TV)"=1/6, "A->G or T->C (TI)"=1/6, "A->T or T->A (TV)"=1/6,
                                   "G->A or C->T (TI)"=1/6, "G->C or C->G (TV)"=1/6, "G->T or C->A (TV)"=1/6),
                               fileType="MGI", dataOut=TRUE, progress=FALSE))
    expect_equal(levels(x[["main"]]$trans_tranv), levels(x[["expect"]]$trans_tranv))
})