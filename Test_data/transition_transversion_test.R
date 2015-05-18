setwd('~/Desktop')

data <- read.delim('master_snvs_indels_somatic.anno.rcnt.filtered.nosilent.tsv')

data <- data[,c('reference', 'variant', 'Extraction_id')]

colnames(data) <- c('reference', 'variant', 'sample')

transition_transversion_plot(data)