test_that("TvTi_convMaf Detects Homozygous events and only counts them once", {
    x <- data.frame(Tumor_Sample_Barcode=c("A"), Reference_Allele="C", Tumor_Seq_Allele1=c("T"), Tumor_Seq_Allele2=c("T"))
    expect <- data.frame(sample=c("A"), reference=c("C"), variant=c("T"))
    expect_equal(TvTi_convMaf(x), expect)
})

test_that("TvTi_convMaf eliminates calls where the reference allele equals the variant allele", {
    x <- data.frame(Tumor_Sample_Barcode=c("A"), Reference_Allele="C", Tumor_Seq_Allele1=c("C"), Tumor_Seq_Allele2=c("C"))
    expect_equal(nrow(TvTi_convMaf(x)), 0)
    
    x <- data.frame(Tumor_Sample_Barcode=c("A"), Reference_Allele="C", Tumor_Seq_Allele1=c("T"), Tumor_Seq_Allele2=c("C"))
    expect_equal(nrow(TvTi_convMaf(x)), 1)
    
    x <- data.frame(Tumor_Sample_Barcode=c("A"), Reference_Allele="C", Tumor_Seq_Allele1=c("C"), Tumor_Seq_Allele2=c("T"))
    expect_equal(nrow(TvTi_convMaf(x)), 1)
})

test_that(TvTi_convMaf)