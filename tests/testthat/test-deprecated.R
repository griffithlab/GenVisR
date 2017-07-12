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

test_that("waterfall_calcMutFreq correctly counts mutations", {
    x <- data.frame(sample=c('a', 'a', 'a', 'b', 'b'), gene=rep('b', 5), trv_type=rep('silent', 5))
    out <- waterfall_calcMutFreq(x)
    expect_equal(out[out$sample == 'a' & out$trv_type == 'Synonymous', 'mutation_total'], 3)
    expect_equal(out[out$sample == 'b' & out$trv_type == 'Synonymous', 'mutation_total'], 2)
    expect_equal(out[out$sample == 'a' & out$trv_type == 'Non Synonymous', 'mutation_total'], 0)
    expect_equal(out[out$sample == 'b' & out$trv_type == 'Non Synonymous', 'mutation_total'], 0)
})

test_that("waterfall_calcMutFreq correctly identifies Synonymous, Non Synonymous mutations", {
    x <- data.frame(sample=c('a', 'a', 'a'), gene=rep('b', 3), trv_type=c('not silent', 'Silent', 'siLenT'))
    out <- waterfall_calcMutFreq(x)
    expect_equal(out[out$sample == 'a' & out$trv_type == 'Synonymous', 'mutation_total'], 2)
    expect_equal(out[out$sample == 'a' & out$trv_type == 'Non Synonymous', 'mutation_total'], 1)   
})

test_that("waterfall_geneAlt checks for proper input", {
    x <- data.frame(sample='a', gene='b', trv_type='missense')
    genes <- c(1) 
    expect_warning(waterfall_geneAlt(x, genes), "not a character vector")
})

test_that("waterfall_geneAlt checks that genes are available for removal", {
    x <- data.frame(sample='a', gene='b', trv_type='missense')
    genes <- c('z')
    expect_warning(waterfall_geneAlt(x, genes), "element not found in x")
})

test_that("waterfall_geneAlt successfully removes genes from input", {
    x <- data.frame(sample=c('a','a'), gene=c('b','c'), trv_type=rep('missense',2))
    genes <- c('b')
    out <- waterfall_geneAlt(x, genes)
    expect_equal(nrow(out), 1)
})

test_that("waterfall_geneRecurCutoff properly detects if recurrence cutoff exceeds max recurrence", {
    x <- data.frame(sample=factor(seq(1:11)), gene=c(rep(c('b', 'c'), 5), 'd'), trv_type=rep('missense', 11))
    recurrence_cutoff <- .51
    expect_warning(waterfall_geneRecurCutoff(x, recurrence_cutoff), "cutoff specified exceeds")
})

test_that("waterfall_geneRecurCutoff properly removes entries not meeting a recurrence cutoff", {
    x <- data.frame(sample=factor(seq(1:11)), gene=c(rep(c('b', 'c'), 5), 'd'), trv_type=rep('missense', 11))
    recurrence_cutoff <- .51
    out <- suppressWarnings(waterfall_geneRecurCutoff(x, recurrence_cutoff))
    expect_equal(nrow(out), 10)    
})

test_that("waterfall_geneSort properly orders genes", {
    x <- data.frame(sample=rep('a', 10), gene=c(rep(c('b', 'c', 'd'), 3), 'e'), trv_type=rep('missense', 10))
    out <- waterfall_geneSort(x) 
    expec <- c('e', 'b', 'c', 'd')
    expect_equal(out, expec)
})

test_that("waterfall_geneSort defers to input to geneOrder if specified by the user", {
    x <- data.frame(sample=rep('a', 10), gene=c(rep(c('b', 'c', 'd'), 3), 'e'), trv_type=rep('missense', 10))
    gene_order <- c("d", "b", "c", "e")
    out <- waterfall_geneSort(x, geneOrder=gene_order)
    expec <- rev(gene_order)
    expect_equal(out, expec)
    
    gene_order <- c("none", "b", "d")
    expec <- c("d", "b")
    out <- suppressWarnings(waterfall_geneSort(x, geneOrder=gene_order))
    expect_equal(out, expec)
})

test_that("waterfall_geneSort warns if it finds genes supplied to geneOrder not in the primary input", {
    # test case where no genes in geneOrder are in primary input
    x <- data.frame(sample=rep('a', 10), gene=c(rep(c('b', 'c', 'd'), 3), 'e'), trv_type=rep('missense', 10))
    gene_order <- c("none")
    expect_warning(waterfall_geneSort(x, geneOrder=gene_order), "Did not find any genes")
    
    # test case where some but not all genes in geneOrder are in primary input
    x <- data.frame(sample=rep('a', 10), gene=c(rep(c('b', 'c', 'd'), 3), 'e'), trv_type=rep('missense', 10))
    gene_order <- c("none", "b", "d")
    expect_warning(waterfall_geneSort(x, geneOrder=gene_order), "The following genes were not found")
})

test_that("waterfall_MAF2anno checks for proper column names",{
    # Check without label column
    x <- data.frame(Tumor_Sample_Barcode='samp1', incorrect='egfr', Variant_Classification='missense')
    label_col <- NULL
    expect_error(waterfall_MAF2anno(x, label_col), "Did not detect")
    
    # Check with label column
    x <- data.frame(Tumor_Sample_Barcode='samp1', Hugo_Symbol='egfr', Variant_Classification='missense', incorrect='R504K')
    label_col <- 'label'
    expect_error(waterfall_MAF2anno(x, label_col), "Did not detect")   
})

test_that("waterfall_MAF2anno output format is correct", {
    # check without label column
    x <- data.frame(Tumor_Sample_Barcode='samp1', Hugo_Symbol='egfr', Variant_Classification='missense')
    label_col <- NULL
    out <- waterfall_MAF2anno(x, label_col)
    expect_equal(sort(colnames(out)), sort(c('sample', 'gene', 'trv_type')))
    
    # check with label column
    x <- data.frame(Tumor_Sample_Barcode='samp1', Hugo_Symbol='egfr', Variant_Classification='missense', label='R504K')
    label_col <- 'label'
    out <- waterfall_MAF2anno(x, label_col)
    expect_equal(sort(colnames(out)), sort(c('sample', 'gene', 'trv_type', label_col)))
    
    # Check class
    expect_is(out, "data.frame")
})

test_that("waterfall_MGI2anno checks for proper column names",{
    # Check without label column
    x <- data.frame(incorrect='samp1', gene_name='egfr', trv_type='missense')
    label_col <- NULL
    expect_error(waterfall_MGI2anno(x, label_col), "Did not detect")
    
    # Check with label column
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense', incorrect='R504K')
    label_col <- 'label'
    expect_error(waterfall_MGI2anno(x, label_col), "Did not detect")   
})

test_that("waterfall_MGI2anno output format is correct", {
    # check without label column
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    label_col <- NULL
    out <- waterfall_MGI2anno(x, label_col)
    expect_equal(sort(colnames(out)), sort(c('sample', 'gene', 'trv_type')))
    
    # check with label column
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense', label='R504K')
    label_col <- 'label'
    out <- waterfall_MGI2anno(x, label_col)
    expect_equal(sort(colnames(out)), sort(c('sample', 'gene', 'trv_type', label_col)))
    
    # Check class
    expect_is(out, "data.frame")
})

test_that("waterfall_NA2gene properly assigns NA values in the input to the top gene", {
    x <- data.frame(sample=c('a', 'b', NA), gene=c('c', 'c', NA), trv_type=c('missense', 'missense', NA))
    out <- waterfall_NA2gene(x)
    expect_false(any(is.na(out$gene)))
})

test_that("waterfall_qual returns a list with three elements", {
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- NULL
    z <- NULL
    file_type <- 'MGI'
    label_col <- NULL
    out <- waterfall_qual(x, y, z, file_type, label_col)
    
    expect_equal(length(out), 3)
    expect_is(out, 'list')
})

test_that("waterfall_qual checks input to file_type is valid", {
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- NULL
    z <- NULL
    file_type <- 'INCORRECT'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col), "Unrecognized")
})

test_that("waterfall_qual verifies input is of proper class", {
    
    # Check x
    x <- c(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- NULL
    z <- NULL
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    # Check y
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- c(sample='samp1', variable='gender', value='male') 
    z <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    # Check z
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- NULL
    z <- c(sample='samp1', mut_burden='5')
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
})

test_that("waterfall_qual verifies correct columns are in input", {
    
    # Check x
    x <- data.frame(INCORRECT='samp1', gene_name='egfr', trv_type='missense',
                    chromosome_name=1, start=1, stop=1, reference="A",
                    variant="T") 
    y <- NULL
    z <- NULL
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    # Check y
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- data.frame(incorrect='samp1', variable='gender', value='male')
    z <- NULL
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- data.frame(sample='samp1', incorrect='gender', value='male')
    z <- NULL
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- data.frame(sample='samp1', variable='gender', incorrect='male')
    z <- NULL
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    # Check z
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- NULL
    z <- c(incorrect='samp1', mut_burden='5')
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- NULL
    z <- c(sample='samp1', incorrect='5')
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')    
})

test_that("waterfall_rmvSilent removes silent mutations", {
    x <- data.frame(sample=rep('a', 3), gene=rep('b', 3), trv_type=c('Silent', 'siLenT', 'missense'))
    out <- droplevels(waterfall_rmvSilent(x))    
    expec <- data.frame(sample=rep('a', 3), gene=c(NA, NA, 'b'), trv_type=c(NA, NA, 'missense'))
    expect_equal(out, expec)
})

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

test_that("waterfall_sampSort properly sorts samples based on a hierarchy", {
    x <- data.frame(sample=c('a', 'a', 'a', 'b', 'b', 'c'), gene=c('x', 'y', 'z', 'x', 'y', 'z'), trv_type=rep('missense', 6))
    out <- waterfall_sampSort(x)
    expec <- c('a', 'c', 'b')
    expect_equal(out, expec)
    
    # Test that convert to boolean properly works on data frame with one gene
    x <- data.frame(sample=c('a', 'b', 'c'), gene=c('x', 'x', 'x'), trv_type=rep('missense', 3))
    out <- waterfall_sampSort(x)
    expec <- c('a', 'b', 'c')
    expect_equal(out, expec)
    
    # Test that samples added in via the plotSamples parameter in waterfall are
    # put last
    x <- data.frame(sample=c('TESTA', 'A'), gene=c(NA, 'x'), trv_type=c(NA, 'missense'))
    out <- waterfall_sampSort(x)
    expec <- c('A', 'TESTA')
    #expect_equal(out, expec)
})

test_that("waterfall_sampSort removes samples in sampOrder not found in primary input", {
    x <- data.frame(sample=c('a', 'a', 'a', 'b', 'b', 'c'), gene=c('x', 'y', 'z', 'x', 'y', 'z'), trv_type=rep('missense', 6))
    sample_order <- c("a", "b", "c", "d")
    out <- suppressWarnings(waterfall_sampSort(x, sampOrder = sample_order))
    expec <- c('a', 'b', 'c')
    expect_equal(out, expec)
})

test_that("waterfall_hierarchyTRV leaves only one entry for a gene/mutation", {
    # Check MGI file
    x <- data.frame(sample=rep(c('a', 'b'), 4), gene=rep(c('c', 'd'), 4), trv_type=rep(c('nonsense', 'missense', 'rna', 'silent'), 2))
    file_type <- 'MGI'
    variant_class_order <- NULL
    out <- waterfall_hierarchyTRV(x, file_type, variant_class_order)
    expect_equal(nrow(out), 2)
    expect_equal(as.character(out[out$sample == 'a' & out$gene == 'c', 'trv_type']), 'nonsense')
    expect_equal(as.character(out[out$sample == 'b' & out$gene == 'd', 'trv_type']), 'missense')
    
    # Check MAF file
    x <- data.frame(sample=rep(c('a', 'b'), 4), gene=rep(c('c', 'd'), 4), trv_type=rep(c('Splice_Site', 'IGR', 'Silent', 'Nonsense_Mutation'), 2))
    file_type <- 'MAF'
    out <- waterfall_hierarchyTRV(x, file_type, variant_class_order)
    expect_equal(nrow(out), 2)
    expect_equal(as.character(out[out$sample == 'a' & out$gene == 'c', 'trv_type']), 'Splice_Site')
    expect_equal(as.character(out[out$sample == 'b' & out$gene == 'd', 'trv_type']), 'Nonsense_Mutation')    
})

test_that("waterfall_hierarchyTRV recognizes an invalid mutation type", {
    # Check for MGI
    x <- data.frame(sample=rep(c('a', 'b'), 4), gene=rep(c('c', 'd'), 4), trv_type=rep(c('missense', 'rna', 'silent', 'invalid'), 2))
    file_type <- 'MGI'
    variant_class_order <- NULL
    expect_error(waterfall_hierarchyTRV(x, file_type, variant_class_order), "invalid mutation type")
    
    # Check for MAF
    x <- data.frame(sample=rep(c('a', 'b'), 4), gene=rep(c('c', 'd'), 4), trv_type=rep(c('Splice_Site', 'IGR', 'Silent', 'invalid'), 2))
    file_type <- 'MAF'
    expect_error(waterfall_hierarchyTRV(x, file_type, variant_class_order), "invalid mutation type")
})