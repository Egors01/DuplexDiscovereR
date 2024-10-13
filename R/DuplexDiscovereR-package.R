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
#' @importFrom GenomicRanges GRanges
#' @import Gviz
#' @import grid
#' @importFrom Biostrings readBStringSet RNAStringSet letterFrequency
#' @importFrom ggsci pal_igv
#' @author Egor Semenchenko
#' @seealso [DuplexDiscovereR vignette](`browseVignettes("DuplexDiscovereR")`)
"_PACKAGE"

# Quiets R CMD CHECK NOTES from tidyverse syntax and func from pkg imported with
# @import
varlist <- c(
    "A", "AB", "A_ovl", "A_span", "B", "B_ovl", "B_span",
    "Pa", "Pb", "ambig", "ambig.A", "ambig.B", "ambig_list.A", "ambig_list.B",
    "arc_max_y", "arm", "artCIGAR", "artCIGAR1", "artCIGAR2", "artSEQ",
    "bad_junction", "brkpt_acceptorB", "brkpt_donorA", "chr_acceptorB",
    "chr_donorA", "chromA", "chromB", "cigar_alnA", "cigar_alnB",
    "cigar_str", ".handleComposite", "cmp", "col11", "col1read1",
    "col1read2", "col2", "col5", "col7", "col8", "col9", "dg_id",
    "dg_id_raw", "distance", "duplex_id",
    "end1", "end2", "endA", "endB", "feature_id",
    "findOverlapPairs", "subsetByOverlaps", "GRanges",
    "mcols<-", "punion", "pintersect", "subjectHits", "IRanges",
    "full_in_plot", "gap", "gene_count.A", "gene_count.B",
    "gi_trans", "grp", "grpvar", "has_A", "has_B", "has_p",
    "hg38_seqlengths", "idcol.1",
    "idcol.2", "in_range", "junction_type", "lengapsA", "lengapsB",
    "map_type", "minus_p", "multigap", "multimap", "nA", "nB",
    "n_i", "n_m", "n_n", "n_p", "n_reads", "n_reads_dg", "n_s",
    "new_duplex_id", "ngapsA", "ngapsB", "nr", "num_chim_aln",
    "old_duplex_id", "overlap_w", "p.adj", "pair_arm",
    "queryHits", "rDens", "ratio.A", "ratio.B", "read_id",
    "read_name", "read_type", "readname", "region_id", "rname",
    "seq_anchor1", "seq_anchor2", "seqlengths", "seqlengths",
    "seqnames1", "seqnames2", "sign_p", "splicejnc", "start1", "start2",
    "startA", "startB", "start_alnA", "start_alnB", "strand1", "strand2",
    "strandA", "strandB", "strand_acceptorB", "strand_donorA",
    "tag", "this_chim_aln_score", "seqlengths<-", ".",
    "vert_id", "vert_id1", "vert_id2", "weight", "width1", "width2",
    "x1", "x2", "x3"
)

if (getRversion() >= "2.15.1") utils::globalVariables(varlist, add = FALSE)
