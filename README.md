# DuplexDiscovereR

<!-- badges: start -->
<!-- badges: end -->

A package for the analysis of the data from RNA duplex probing experiments
## About

DuplexDiscovereR is a package for analysing data from RNA cross-linking and 
proximity ligation protocols such as SPLASH, PARIS, LIGR-seq and other
experimental methods that capture information about RNA-RNA interactions 
as chimeric fragments after high-throughput sequencing.
It provides functions for uniform processing and visualisation of the raw chimeric 
read data. The methods implemented in DuplexDiscovereR are compatible with 
Bioconductor's GenomicRanges and InteractionSet classes, allowing for the 
integration of RNA-RNA interaction probing results into transcriptomic
analyses. 

## License

The `DuplexDiscovereR` code is distributed under [GNU Public License 3.0](https://www.gnu.org/licenses/gpl-3.0.html). 
The documentation, including the manual pages and vignettes,
is distributed under a [CC BY-SA 4.0 license](https://creativecommons.org/licenses/by-sa/4.0/deed.en).


Some functionality of `DuplexDiscovereR` depends on the [RNAduplex](https://www.tbi.univie.ac.at/RNA/RNAduplex.1.html) from the
[ViennaRNA](https://www.tbi.univie.ac.at/RNA/). ViennaRNA is distributed under its own licence.


## Installation

To install this package, use R (version "4.4"):

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DuplexDiscovereR")
```

You can also install the package from Github with 
```r
library(devtools)
install_github('Egors01/DuplexDiscovereR')
```
 [RNAduplex](https://www.tbi.univie.ac.at/RNA/RNAduplex.1.html) can be installed by foolowing instructions on [ViennaRNA web-page](https://www.tbi.univie.ac.at/RNA/)

## Getting started

Package tutorial is avaliable on [Bioconductor](https://www.bioconductor.org/packages/release/bioc/vignettes/DuplexDiscovereR/inst/doc/DuplexDiscovereR.html) 

It can be loaded from R by 

```r
browseVigenttes("DuplexDiscovereR")
```


## Citataion
Semenchenko E, Tsybulsky V, Meyer IM (2025). "DuplexDiscovereR: a computational method for the analysis of experimental duplex RNA-RNA interaction data". Nucleic Acids Research. 
https://doi.org/10.1093/nar/gkaf266
```
@article{DuplexDiscovereR,
  title={DuplexDiscovereR: a computational method for the analysis of experimental duplex RNA-RNA interaction data},
  author={Semenchenko E, Tsybulsky V, Meyer IM},
  journal={Nucleic Acids Research},
  doi = {10.1093/nar/gkaf266}
  year={2018},
  publisher={Oxford University Press}
}
```


