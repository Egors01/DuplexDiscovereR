#' DuplexDiscoverer
#' DuplexDiscoverer is a package for analysis of the data from RNA cross-linking and proximity
#' ligation protocols such as SPLASH, PARIS, LIGR-seq and others.
#' DuplexDiscoverer takes input in a form of Chimericly or split -aligned reads.
#' It implements procedures for classification, filtering and efficient clustering
#' of inidividual read alignments into duplex groups (DGs)
#' After DGs are found, RNA duplex formation and their hybridization energies
#' are predicted,additional metrics
#' i.e p-values or mean alignement scores can be calcualted to rank and analyze
#' final RNA duplex predictions.Data from multiple experiments or replicates
#' can be processed separately and further compared to analyse reproducability
#' of experimental method.
#'
#' @name DuplexDiscoverer
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
#' @seealso [`browseVignettes("DuplexDiscoverer")`
"_PACKAGE"
