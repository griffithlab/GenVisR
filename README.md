[![Build Status](https://travis-ci.org/griffithlab/GenVisR.svg?branch=master)](https://travis-ci.org/griffithlab/GenVisR)

## GenVisR

Please cite: "Skidmore et al. 2016 GenVisR: Genomic Visualizations in R Bioinformatics 32, 3012-3014" [pubmed](https://www.ncbi.nlm.nih.gov/pubmed/27288499)

[Bioconductor](https://bioconductor.org/packages/release/bioc/html/GenVisR.html)

Intuitively visualizing and interpreting data from high-throughput genomic technologies continues to be challenging. "Genomic Visualizations in R" (GenVisR) attempts to alleviate this burden by providing highly customizable publication-quality graphics supporting multiple species and focused primarily on a cohort level (i.e., multiple samples/patients). GenVisR attempts to maintain a high degree of flexibility while leveraging the abilities of ggplot2 and bioconductor to achieve this goal.

### Install from Bioconductor

For the majority of users we recommend installing GenVisR from the release branch of Bioconductor, Installation instructions using this method can be found on the [GenVisR landing page](http://bioconductor.org/packages/GenVisR/) on Bioconductor.

Please note that GenVisR imports a few packages that have "system requirements", in most cases these requirements will already be installed. If they are not please follow the instructions to install these packages given in the R terminal. Briefly these packages are: "libcurl4-openssl-dev" and "libxml2-dev"

### Development

Development for GenVisR occurs on the griffith lab github repository available [here](https://github.com/griffithlab/GenVisR). For users wishing to contribute to development we recommend cloning the GenVisR repo there and submitting a pull request. Please note that development occurs on the R version that will be available at each Bioconductor release cycle. This ensures that GenVisR will be stable for each Bioconductor release but it may necessitate developers download R-devel.

We also encourage users to report bugs and suggest enhancements to GenVisR on the github issue page available [here](https://github.com/griffithlab/GenVisR/issues):

To install the latest development version of GenVisR (not recommended for most users):

```R
# install and load devtools package
install.packages("devtools")
library(devtools)

# install GenVisR from github
install_github("griffithlab/GenVisR")
```

### Vignettes

Documentation for GenVisR can be found on the bioconductor landing page in the form of vignettes available here [GenVisR](http://bioconductor.org/packages/GenVisR/). Tutorials can also be found on biostars.org. Vignettes can also be viewed from within R.

```R
# view vignettes
vignette(package="GenVisR")
```
