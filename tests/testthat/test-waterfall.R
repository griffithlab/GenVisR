test_that("waterfall levels data correctly and as such data aligns", {
    out <- waterfall(brcaMAF, plotGenes=c("PIK3CA", "TP53"), out="data")
    expect_equal(levels(out[["main"]]$sample), levels(out[["mutation_count"]]$sample))    
})

test_that("waterfall correctly counts mutations in a sample", {
    # Obtain output
    out <- waterfall(brcaMAF, plotGenes=c("PIK3CA", "TP53"), out="data")
    out_A <- out$mutation[out$mutation_count$sample == "TCGA-A1-A0SI-01A-11D-A142-09" &
                              out$mutation_count$trv_type == "Synonymous",]
    out_B <- out$mutation[out$mutation_count$sample == "TCGA-A1-A0SI-01A-11D-A142-09" &
                              out$mutation_count$trv_type == "Non Synonymous",]
    
    # Calculate expected output
    expec <- brcaMAF[brcaMAF$Tumor_Sample_Barcode=="TCGA-A1-A0SI-01A-11D-A142-09",]
    expec_A <- nrow(expec[expec$Variant_Classification == "Silent" ,c(1, 9, 16)])
    expec_B <- nrow(expec[expec$Variant_Classification != "Silent" ,c(1, 9, 16)])
    
    expect_equal(out_A$mutation_total, expec_A)
    expect_equal(out_B$mutation_total, expec_B)
})