# DuplexDiscoverer

<!-- badges: start -->
<!-- badges: end -->

A package for the analysis of the data from RNA duplex probing experiments
## About

DuplexDiscoverer is a package for analysing data from RNA cross-linking and 
proximity ligation protocols such as SPLASH, PARIS, LIGR-seq and other
experimental methods that record information about RNA-RNA interactions 
as chimeric fragments after high-throughput sequencing.
It provides functions for unified processing and visualisation of the raw chimeric 
read data. The methods implemented in DuplexDiscoverer are compatible with 
Bioconductor's GenomicRanges and InteractionSet classes, allowing the 
integration of RNA-RNA interaction probing results into transcriptomic
analysis. 

## License

The `DuplexDiscoverer` code is distributed under [GNU Public License 3.0](https://www.gnu.org/licenses/gpl-3.0.html). 

Some functionality of `DuplexDiscoverer` depends on the [RNAduplex](https://www.tbi.univie.ac.at/RNA/RNAduplex.1.html) from the
[ViennaRNA](https://www.tbi.univie.ac.at/RNA/). ViennaRNA is distributed under its own licence.


## Installation

To install this package, start R (version "4.4") and enter:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DuplexDiscoverer")
```

You can install the development version from Github with 
```r
library(devtools)
install_github('Egors01/DuplexDiscoverer')
```
 [RNAduplex](https://www.tbi.univie.ac.at/RNA/RNAduplex.1.html) can be installed by foolowing instructions on [ViennaRNA web-page](https://www.tbi.univie.ac.at/RNA/)

## Getting started

See the tutorial by loading package vignette 
```r
browseVigenttes("DuplexDiscoverer")
```
## Citataion

*Not avaliable before paper publication*

