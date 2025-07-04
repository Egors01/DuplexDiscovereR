# Version  0.99.0 
- Date: 2024-06-09
- Submission to Bioconductor

# Version  0.99.1
- Date: 2024-10-11
- Modifications according to Bioconductor - review requirements
- Fixed bugs in deduplication function
- Changed splice junction detection getSpliceJunctionChimeras() to account for 
alt 3' and 5' ends

# Version  1.1.1
- Date: 2024-12-30
- Added added raw abundance values and duplex counts for 2x2 contingency table to the p-value calculation output

# Version  1.1.2
- Date: 2025-03-12
- Removed BiocCheck files 

# Version  1.2.1
- Date: 2025-07-04
- Added trimming of the PE alignments relative to the ligation point. 
    - For the PE reads processed by STAR, the type of the junction is detected. 
    - To narrow down the hybrid prediction, user can trim reads - N_trim nt will be left adjacent to the ligation point.
- Added missing arguments to the function `runDuplexDiscoverer` which runs the whole pipeline. Keys for annotation, trimming and minimum chimeric length cutoffs can be set by user and will be passed to the functions called from inside.

