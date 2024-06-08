#' Chimeric reads of SPLASH
#'
#' A Chimeric.out.Junction file with a subset of chr 22 Chimeric reads
#' detected by SPLASH protocol in Human embryonic stem cells.
#' @source [SequenceReadArcive](https://www.ncbi.nlm.nih.gov/sra/?term=SRR3404943)
#' Reads were aligned with STAR
#' see \code{system.file("extdata/scripts", "DD_data_generation.R", package =
#' "DuplexDiscovereR")} for details on the pre-processing and sub-setting the
#' data
#' @docType data
#' @keywords datasets
#' @usage data(RNADuplexesSampleData)
#' @returns `tibble` with columns of Chimeric.junction.out
"RNADuplexesRawChimSTAR"

#' Chimeric reads of SPLASH converted to .bedpe fromat
#'
#' A Chimeric.out.Junction file with a subset of chr 22 Chimeric reads
#' detected by SPLASH protocol in Human embryonic stem cells.
#' @source [SequenceReadArcive](https://www.ncbi.nlm.nih.gov/sra/?term=SRR3404943)
#' Reads were aligned with STAR and filtered to contain only reads
#' which could be represented as 2-arm chimeric alignments. Converted to the
#' [bedpe](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format)
#' format
#' see \code{system.file("extdata/scripts", "DD_data_generation.R", package =
#' "DuplexDiscovereR")} for details on the pre-processing and sub-setting the
#' data
#' @docType data
#' @keywords datasets
#' @usage data(RNADuplexesSampleData)
#' @returns `tibble` with columns of bedpe format
"RNADuplexesRawBed"

#' RNA duplex reads of SPLASH, clustered and assigned to duplex groups
#'
#' `GInteractions` read-level object containing processed reads,annotated with duplex group
#'  ids, read types gene names and p-values
#' @source [SequenceReadArcive](https://www.ncbi.nlm.nih.gov/sra/?term=SRR3404943)
#' Reads were aligned with STAR and duplex groups were identified
#' see \code{system.file("extdata/scripts", "DD_data_generation.R", package =
#' "DuplexDiscovereR")} for details on the data generation proccedure.
#' @docType data
#' @keywords datasets
#' @usage data(RNADuplexesSampleData)
#' @returns `GInteractions` with
#' - `n_reads_dg ` : number of reads in the duplex group (DG)
#' - `duplex_id  ` : temporary id for RNA duplexes which could be found before
#' clustering (duplicated or shifted by couple of nt )
#' - `dg_id` :id of the duplex group
#' - `score` : median alignment score in duplex group
#' - other columns inherited from the STAR Chimeric.out.Junction
"RNADuplexSampleClustReads"


#' RNA duplex reads of SPLASH derived from chimeric alignments
#'
#' `GInteractions` read-level object containing two-arm chimeric reads extracted from
#' mapping output and which can be represented in the `GInteraction` object
#' @source [SequenceReadArcive](https://www.ncbi.nlm.nih.gov/sra/?term=SRR3404943)
#' @details
#' see \code{system.file("extdata/scripts", "DD_data_generation.R", package =
#' "DuplexDiscovereR")} for details on the data generation proccedure.
#' @docType data
#' @keywords datasets
#' @usage data(RNADuplexesSampleData)
#' @returns `GInteractions` with
#' - `readname` : read name
#' - `map_type ` : type of the mapped read (2arm by design of pre-filtering)
#' - `junction_type` : if read jucntion is too short, or it not a 'true' ligated
#' reads because of the jucntoin coincides with splice junction
#' - cigar_aln* columns inherited from the STAR Chimeric.out.Junction output
"RNADuplexSampleGI"



#' RNA duplex reads of SPLASH, clustered and collapsed to duplex groups
#'
#' `GInteractions` duplex group -level object containing detected duplex groups,
#' annotated with duplex group ids, gene_names and p-values
#' @source [SequenceReadArcive](https://www.ncbi.nlm.nih.gov/sra/?term=SRR3404943)
#' Reads were aligned with STAR and duplex groups were identified
#' see \code{system.file("extdata/scripts", "DD_data_generation.R", package =
#' "DuplexDiscovereR")} for details on the data generation procedure.
#' @docType data
#' @keywords datasets
#' @usage data(RNADuplexesSampleData)
#' @returns `GInteractions` with
#' - `n_reads` : number of reads in the duplex group (DG)
#' - `dg_id` :id of the duplex group
#' - `p_val` : BH adjusted p-value of testing to reject hypothesis of DG arising
#'  from random ligation
#' - `score` : median alignment score in duplex group
#' - other columns with `.A` and `.B` annotating to which genes either arm of the
#' DG maps
"RNADuplexSampleDGs"


#' Gene coordinates on human chromosome 22
#'
#' `Granges` containing gene coordinates of human chromosome 22 obtained from
#' GENCODEv44 annotaion
#' @source [GENCODEv44](https://www.gencodegenes.org/human/release_44.html)
#' @details
#' see \code{system.file("extdata/scripts", "DD_data_generation.R", package =
#' "DuplexDiscovereR")} for details
#' @docType data
#' @keywords datasets
#' @usage data(RNADuplexesSampleData)
#' @returns statdatd GENCODE gtf fields
"SampleGeneAnnoGR"


#' Gene coordinates on human chromosome 22
#'
#' `Granges` containing coordinates of splice junctions human chromosome 22 obtained from
#' GENCODEv44 annotaion
#' @source [GENCODEv44](https://www.gencodegenes.org/human/release_44.html)
#' @details
#' see \code{system.file("extdata/scripts", "DD_data_generation.R", package =
#' "DuplexDiscovereR")} for details
#' @docType data
#' @keywords datasets
#' @usage data(RNADuplexesSampleData)
#' @returns statdatd GENCODE gtf fields
"SampleSpliceJncGR"


#'  Gene counts on human chromosome 22, embryonic stem cells
#'
#' File generated by mapping with STAR using `--quantMode GeneCounts`
#' see \code{system.file("extdata/scripts", "DD_data_generation.R", package =
#' "DuplexDiscovereR")} for details on the pre-processing and sub-setting the
#' @source [SequenceReadArcive](https://www.ncbi.nlm.nih.gov/sra/?term=SRR3404943)
#' @docType data
#' @keywords datasets
#' @usage data(RNADuplexesSampleData)
#' @returns `tibble` with columns of Chimeric.junction.out
"RNADuplexesGeneCounts"

#' RNA duplex reads of SPLASH derived from chimeric alignments
#'
#' `GInteractions` object containing two-arm chimeric reads extracted from
#' mapping output and which can be represented in the `GInteraction` object
#' and subset to chr22: 23877144-45562960 '*'
#' @source [SequenceReadArcive](https://www.ncbi.nlm.nih.gov/sra/?term=SRR3404943)
#' @details
#' see \code{system.file("extdata/scripts", "DD_data_generation.R", package =
#' "DuplexDiscovereR")} for details on the data generation procedure.
#' @docType data
#' @keywords datasets
#' @usage data(RNADuplexesSmallGI)
#' @name SampleSmallGI
"SampleSmallGI"
