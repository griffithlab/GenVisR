# GGgenome
Intuitively visualizing and interpreting data from high-throughput genomic technologies continues to be challenging. GGgenome attempts to alleviate this burden by providing highly customizable publication-quality graphics focused primarily on a cohort level (i.e., multiple samples/patients). GGgenome attempts to maintain a high degree of flexibility while leveraging the abilities of ggplot2 and bioconductor to achieve this goal.

## Installation Instructions

* Install devtools
```R
install.packages("devtools")
```

* Install Bioconductor dependencies
```R
source("http://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges", "biomaRt", "GenomicFeatures", "UniProt.ws"))
```

* It is suggested but not required to install knitr and rmarkdown for vignette creation
```R
install.packages(c("rmarkdown", "knitr"))
```

* Install GGgenome
```R
devtools::install_github("griffithlab/Gggenome")
```
