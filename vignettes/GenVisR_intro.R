## ----kable 1, echo=FALSE, error=TRUE-------------------------------------
library(knitr)
MGI <- c("nonsense", "frame_shift_del",
         "frame_shift_ins", "splice_site_del",
         "splice_site_ins", "splice_site",
         "nonstop", "in_frame_del", "in_frame_ins",
         "missense", "splice_region_del",
         "splice_region_ins", "splice_region",
         "5_prime_flanking_region",
         "3_prime_flanking_region",
         "3_prime_untranslated_region",
         "5_prime_untranslated_region", "rna",
         "intronic", "silent")
MAF <- c("Nonsense_Mutation", "Frame_Shift_Ins",
         "Frame_Shift_Del", "Translation_Start_Site",
         "Splice_Site", "Nonstop_Mutation",
         "In_Frame_Ins", "In_Frame_Del",
         "Missense_Mutation", "5\'Flank",
         "3\'Flank", "5\'UTR", "3\'UTR", "RNA", "Intron",
         "IGR", "Silent", "Targeted_Region", "", "")

kable(as.data.frame(cbind(MAF, MGI)))

## ---- eval=FALSE, error=TRUE---------------------------------------------
#  # Plot the mutation landscape
#  waterfall(brcaMAF, fileType="MAF")

## ---- fig.keep='last', fig.width=10, fig.height=7, message=FALSE, warning=FALSE, results='hide', tidy=TRUE, error=TRUE----
# Load GenVisR and set seed
library(GenVisR)
set.seed(383)

# Plot only genes with mutations in 6% or more of samples
waterfall(brcaMAF, mainRecurCutoff=.06)

## ---- fig.keep='last', fig.width=10, fig.height=7, message=FALSE, warning=FALSE, results='hide', tidy=TRUE, error=TRUE----
# Plot only the specified genes
waterfall(brcaMAF, plotGenes=c("PIK3CA", "TP53", "USH2A", "MLL3", "BRCA1"))

## ----kable, echo=FALSE, tidy=TRUE, error=TRUE----------------------------
kable(as.data.frame(cbind(sample=as.character(brcaMAF[1:10,16]),
                          mut_burden=as.numeric(rnorm(10, mean=2, sd=.5)))))

## ---- fig.keep='last', fig.width=12, fig.height=8.5, message=FALSE, warning=FALSE, results='hide', tidy=TRUE, error=TRUE----
# Create clinical data
subtype <- c('lumA', 'lumB', 'her2', 'basal', 'normal')
subtype <- sample(subtype, 50, replace=TRUE)
age <- c('20-30', '31-50', '51-60', '61+')
age <- sample(age, 50, replace=TRUE)
sample <- as.character(unique(brcaMAF$Tumor_Sample_Barcode))
clinical <- as.data.frame(cbind(sample, subtype, age))

# Melt the clinical data into "long" format.
library(reshape2)
clinical <- melt(clinical, id.vars=c('sample'))

# Run waterfall
waterfall(brcaMAF, clinDat=clinical,
          clinVarCol=c('lumA'='blue4', 'lumB'='deepskyblue', 
                            'her2'='hotpink2', 'basal'='firebrick2',
                            'normal'='green4', '20-30'='#ddd1e7',
                            '31-50'='#bba3d0', '51-60'='#9975b9',
                            '61+'='#7647a2'), 
          plotGenes=c("PIK3CA", "TP53", "USH2A", "MLL3", "BRCA1"),
          clinLegCol=2,
          clinVarOrder=c('lumA', 'lumB', 'her2', 'basal', 'normal',
                         '20-30', '31-50', '51-60', '61+'))

