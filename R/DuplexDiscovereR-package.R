#' DuplexDiscovereR
#' 
#' @title Analysis of the data from RNA duplex probing experiments
#' @description
#' DuplexDiscovereR is a package for analysing data from RNA cross-linking and 
#' proximity ligation protocols such as SPLASH, PARIS, LIGR-seq and others, 
#' which provide information about intra-molecular RNA-RNA interactions through 
#' chimeric RNA-seq reads. Chimerically aligned fragments in these experiments 
#' correspond to the base-paired stretches (RNA duplexes) of RNA molecules .
#' DuplexDiscovereR takes input in the form of chimericly or split -aligned reads,
#' It implements procedures for alignment classification, filtering and efficient
#' clustering  of  individual chimeric reads into duplex groups (DGs).
#' Once DGs are found, RNA duplex formation and their hybridization energies are 
#' predicted. Additional metrics, such as p-values or mean DG alignment scores, 
#' can be calculated to rank and analyse the final set of RNA duplexes. 
#' Data from multiple experiments or replicates can be processed separately 
#' and further compared to check the reproducibility of the experimental method. 
#'
#' @name DuplexDiscovereR
#'
#' @import rtracklayer
#' @import InteractionSet
#' @import digest
#' @import rlang
#' @importFrom scales rescale
#' @importFrom GenomicAlignments cigarWidthAlongReferenceSpace
#' @import tibble
#' @import dplyr
#' @import stringr
#' @importFrom tidyr unite replace_na
#' @importFrom igraph graph_from_data_frame simplify decompose cluster_louvain cluster_fast_greedy membership edge.attributes V
#' @importFrom grDevices col2rgb  hsv  rgb2hsv
#' @importFrom methods as callNextMethod getClassDef new is getClass slot
#' @importFrom stats median p.adjust pbinom setNames
#' @importFrom utils head setTxtProgressBar txtProgressBar write.table
#' @importFrom purrr set_names
#' @import Gviz
#' @import grid
#' @importFrom Biostrings readBStringSet RNAStringSet letterFrequency
#' @importFrom ggsci pal_igv
#' @author Egor Semenchenko
#' @seealso [`browseVignettes("DuplexDiscovereR")`
"_PACKAGE"
