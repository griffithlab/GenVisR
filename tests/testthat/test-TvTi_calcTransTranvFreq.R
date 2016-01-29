test_that("TvTi_calcTransTranvFreq adds in trans/tranv missing from samples", {
    x <- data.frame(reference=c("A", "G"),
                    variant=c("T", "C"),
                    sample=c("samp1", "samp2"),
                    trans_tranv=c("A->T or T->A (TV)", "A->G or T->C (TI)"))
    expect_equal(nrow(TvTi_calcTransTranvFreq(x)), 12)
})

test_that("TvTi_calcTransTranvFreq correctly calculates the proportions of TR/TV", {
    x <- data.frame(reference=c("A", "G"),
                    variant=c("T", "C"),
                    sample=c("samp1", "samp1"),
                    trans_tranv=c("A->T or T->A (TV)", "A->G or T->C (TI)")) 
    out <- TvTi_calcTransTranvFreq(x)
    expect_equal(out[out$trans_tranv == "A->G or T->C (TI)",]$Prop, .5)
    
    x <- data.frame(reference=c("A", "G"),
                    variant=c("T", "C"),
                    sample=c("samp1", "samp2"),
                    trans_tranv=c("A->T or T->A (TV)", "A->G or T->C (TI)"))
    out <- TvTi_calcTransTranvFreq(x)
    expect_equal(sum(out$Prop), 2)
})

test_that("TvTi_calcTransTranvFreq correctly calculates the frequency of Tv/Ti", {
    x <- data.frame(reference=c("A", "G"),
                    variant=c("T", "C"),
                    sample=c("samp1", "samp1"),
                    trans_tranv=c("A->T or T->A (TV)", "A->G or T->C (TI)")) 
    out <- TvTi_calcTransTranvFreq(x)
    expect_equal(out[out$trans_tranv == "A->G or T->C (TI)",]$Freq, 1)
    expect_equal(out[out$trans_tranv == "A->T or T->A (TV)",]$Freq, 1)
    
    x <- data.frame(reference=c("A", "G"),
                    variant=c("T", "C"),
                    sample=c("samp1", "samp2"),
                    trans_tranv=c("A->T or T->A (TV)", "A->G or T->C (TI)"))
    out <- TvTi_calcTransTranvFreq(x)
    expect_equal(sum(out$Freq), 2)
})