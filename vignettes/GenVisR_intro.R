## ---- message=FALSE------------------------------------------------------
library(GenVisR)

## ---- fig.keep='last', fig.width=10, fig.height=6.5, message=FALSE, warning=FALSE, results='hide'----
mutSpec(brcaMAF)

## ---- fig.keep='last', fig.width=10, fig.height=6.5, message=FALSE, warning=FALSE, results='hide'----
mutSpec(brcaMAF, main.recurrence_cutoff=3)

## ---- fig.keep='last', fig.width=10, fig.height=6.5, message=FALSE, warning=FALSE, results='hide'----
mutSpec(brcaMAF, main.genes=c("PIK3CA", "TP53", "USH2A", "MLL3", "BRCA1"))

## ----kable, echo=FALSE---------------------------------------------------
library(knitr)
kable(as.data.frame(cbind(sample=as.character(brcaMAF[1:10,16]), mut_burden=as.numeric(rnorm(10, mean=2, sd=.5)))))

## ---- fig.keep='last', fig.width=10.5, fig.height=7.5, message=FALSE, warning=FALSE, results='hide'----
# create fake clinical data
subtype <- c('lumA', 'lumB', 'her2', 'basal', 'normal')
subtype <- sample(subtype, 50, replace=TRUE)
age <- c('20-30', '31-50', '51-60', '61+')
age <- sample(age, 50, replace=TRUE)
sample <- as.character(unique(brcaMAF$Tumor_Sample_Barcode))
clinical <- as.data.frame(cbind(sample, subtype, age))

# melt the data
library(reshape2)
clinical <- melt(clinical, id.vars=c('sample'))

# Run mutSpec
mutSpec(brcaMAF, clinDat=clinical, clin.var.colour=c('lumA'='blue4', 'lumB'='deepskyblue',
'her2'='hotpink2', 'basal'='firebrick2', 'normal'='green4', '20-30'='#ddd1e7',
'31-50'='#bba3d0', '51-60'='#9975b9', '61+'='#7647a2'),
main.genes=c("PIK3CA", "TP53", "USH2A", "MLL3", "BRCA1"), clin.legend.col=2, clin.var.order=
c('lumA', 'lumB', 'her2', 'basal', 'normal', '20-30', '31-50', '51-60', '61+'))

