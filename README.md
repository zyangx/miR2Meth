
# miR2Meth
##Cell Type Deconvolution using miRNA Expression Profliling"

================

<!-- badges: start -->
<!-- badges: end -->

The **miR2Meth** package provides a tool to construct DNA methylation profile for miRNAs by mapping the CpG sites from the widely-used Illumina methylation BeadArray platforms to miRNA promoters. Users can upload the DNA methylatiion matrix data generated by the platform of 450K, EPIC or EPICv2 assay, the miRNA promoter methylation profile is then constructed by using the methylation levels of the probes mapping into genomic regions of miRNA promoters. miR2Meth also provides functions to perform the differential methylation analysis of miRNAs. In addition, by integrating miRNA expression profiling data, users are also able to identify aberrant miRNA that associated with promoter methylation dysregulation. 

## Installation

You can install the development version of miR2Meth like so:

``` r
devtools::install_github("zyangx/miR2Meth")
```

## Example

This is a basic example which shows you how to construct DNA methylation profile for miRNAs:

``` r
library(miR2Meth)
## basic example code

data(beta450k.m)

# obtain methylation profile for miRNAs by a single function
miRBeta.m <- probe2miR(beta=beta450k.m, array="450k", 
                       type="Mature", upstream=2000,
                       downstream=2000, method="mean")

# miRNA methylation matrix with miRNAs as rows and samples as columns 
head(miRBeta.m)
```
