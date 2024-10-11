#' Wrapper for classification of the 2arm chimeric reads
#'
#' @description
#' Wraps two procedures for different types of classification for read alignment:
#'  \describe{
#'   \item{overlap type}{test if chimeric junction map to two non-overlapped regions or shorter than defined minimum distance}
#'   \item{splice junction}{test if chimeric junction is also a splice junction}
#' }
#' @details
#' Calls detection of the chimeric junction type, annotates short junctions
#' on same chromosome an strand as 'short'. Compares chimeric
#' junctions with splice junctions. Adds results as the new metadata fields
#' parallel to the input.
#'
#' @param gi `GInteractions` object
#' @param min_junction_len minimum allowed distance between two chimeric arms
#' @param junctions_gr `Granges` object with the splice junctions coordinates
#' @param max_sj_shift maximum shift between either donor and acceptor splice
#' sites and corresponding chimreic junction coordinates to count chimeric
#' junction as splice junction
#' @returns `GInteractions` object object of the same size with new columns:
#' \describe{
#' \item{splicejnc}{filled with 0 or 1}
#' \item{junction_type}{factor for the junction types}
#' }
#' @export
#' @seealso [DuplexDiscovereR::getChimericJunctionTypes(),
#' DuplexDiscovereR::getSpliceJunctionChimeras()]
#' @examples
#' data("RNADuplexesSampleData")
#' head(RNADuplexSampleGI)
#' # remove all metadata
#' mcols(RNADuplexSampleGI) <- NULL
#' gi <- classifyTwoArmChimeras(RNADuplexSampleGI,
#'     min_junction_len = 5,
#'     junctions_gr = SampleSpliceJncGR, max_sj_shift = 10
#' )
#' table(gi$splicejnc)
#' table(gi$junction_type)
classifyTwoArmChimeras <- function(gi, min_junction_len = 4,
    junctions_gr, max_sj_shift = 4) {
    gi <- getChimericJunctionTypes(gi, normal_gap_threshold = min_junction_len)
    gi <- getSpliceJunctionChimeras(gi, sj_gr = junctions_gr, sj_tolerance = max_sj_shift)

    return(gi)
}


#' Classify chimeric junctions of two-arm reads into types
#' @description
#' Chimeric reads which can be represented ans two-arm interactions can be
#' divided into several categories based on the distance between the chimeric
#' fragments and existence of the overlap between these fragments.
#' @details
#' Takes `GInteractions` object and classifies junctions into following categories
#' \describe{
#'   \item{2arm}{ normal chimeric read }
#'   \item{2arm_short}{ normal chimeric read with junction < \emph{normal_gap_threshold}}
#'   \item{self_ovl}{ arms overlap}
#'   \item{antisense_ovl}{ arms overlap on the opposite strand}
#' }
#' @param gi \code{GInteractions} object
#' @param normal_gap_threshold minimum allowed distance between chimeric arms
#
#' @return gi object of the same size with the 'junction_type' field added
#' @export
#' @examples
#' data("RNADuplexesSampleData")
#' preproc_df <- runDuplexDiscoPreproc(RNADuplexesRawBed, table_type = "bedpe")
#' preproc_gi <- makeGiFromDf(preproc_df)
#' preproc_gi <- getChimericJunctionTypes(preproc_gi)
#' table(preproc_gi$junction_type)
getChimericJunctionTypes <- function(gi, normal_gap_threshold = 10) {
    message("\n--- filtering 2-arm duplexes by junction type  ---")
    gi$junction_type <- factor(x = "2arm", levels = c(
        "2arm",
        "2arm_short",
        "self_ovl",
        "antisense_ovl"
    ))
    gi$idx <- seq_len(length(gi))
    gi$gap <- pairdist(gi, type = "gap")
    gi$intra <- pairdist(gi, type = "intra")
    gi$same_strand <- as.integer(gi@regions[gi@anchor1]@strand == gi@regions[gi@anchor2]@strand)

    expr <- ((gi$intra == TRUE) &
        (gi$same_strand == 1) &
        !is.na(gi$gap) &
        (gi$gap < normal_gap_threshold))
    if (any(expr)) {
        gi[expr]$junction_type <- "2arm_short"
    }

    expr <- ((gi$intra == TRUE) &
        (gi$same_strand == 0) &
        (gi$gap < 0))
    if (any(expr)) {
        gi[expr]$junction_type <- "antisense_ovl"
    }


    expr <- ((gi$intra == TRUE) &
        (gi$same_strand == 1) &
        (gi$gap < 0))
    if (any(expr)) {
        gi[expr]$junction_type <- "self_ovl"
    }

    cts <- table(mcols(gi)$junction_type)

    gi$idx <- NULL
    gi$gap <- NULL
    gi$same_strand <- NULL
    gi$intra <- NULL

    nshort <- ifelse(!is.na(cts["2arm_short"]), cts["2arm_short"], 0)
    nantisense <- ifelse(!is.na(cts["antisense_ovl"]), cts["antisense_ovl"], 0)
    nself_ovl <- ifelse(!is.na(cts["self_ovl"]), cts["self_ovl"], 0)

    message("Duplexes categorized by chimeric junction span: ", length(gi))
    message(
        "Duplexes with normal junctions: ", cts["2arm"], ": ",
        round(cts["2arm"] / length(gi) * 100, 2), " %"
    )
    message(
        "Duplexes with self-overlap: ", nself_ovl, ": ",
        round(nself_ovl / length(gi) * 100, 2), " %"
    )
    message(
        "Duplexes with antisense self-overlap: ", nantisense,
        ": ", round(nantisense / length(gi) * 100, 2), " %"
    )
    message(
        "Duplexes with the junction shorter than ", normal_gap_threshold,
        "nt): ", nshort, ": ",
        round(nshort / length(gi) * 100, 2), " %"
    )

    return(gi)
}


#' Identify chimeric junctions coinciding with the splice junctions
#'
#' @description
#' Marks interactions which starts/ends within specified shift from the
#' known splice junctions.
#' @param gi \pkg{GInteractions} object
#' @param sj_gr \pkg{Granges} object with the splice junctions data
#' @param sj_tolerance maximum shift between either donor and acceptor splice
#' sites and corresponding chimreic junction coordinates to count chimeric
#' junction as splice junction
#'
#' @return gi object with added 'splicejnc' field
#' @export
#' @examples
#' data("RNADuplexesSampleData")
#' gi <- getSpliceJunctionChimeras(RNADuplexSampleGI, SampleSpliceJncGR)
#' table(gi$splicejnc)
getSpliceJunctionChimeras <- function(gi, sj_gr, sj_tolerance = 20) {
    message("\n--- searching for the exon-exon junctions  ---")
    gi$idx <- seq_len(length(gi))
    gi$gap <- pairdist(gi, type = "gap")
    gi$intra <- pairdist(gi, type = "intra")
    gi$same_strand <- as.integer(gi@regions[gi@anchor1]@strand == gi@regions[gi@anchor2]@strand)
    slct <- tibble(x1 = gi$intra == TRUE, x2 = gi$same_strand == 1, x3 = gi$gap >= 1) %>%
        replace(is.na(.), 0) %>%
        mutate(slct = x1 & x2 & x3) %>%
        pull(slct)

    intragi <- gi[slct]
    intergi <- gi[!slct]

    gr_chim <- get_chimeric_junctions_onestrand(intragi)
    gr_chim$idx <- intragi$idx

    SjHits <- findOverlaps(gr_chim, sj_gr, type = "equal", maxgap = sj_tolerance)
    duplexes_coincide_sj <- unique(gr_chim[queryHits(SjHits)]$idx)

    gi$splicejnc <- left_join(tibble("idx" = gi$idx),
        tibble("idx" = duplexes_coincide_sj, "splicejnc" = 1),
        by = "idx"
    ) %>%
        replace(is.na(.), 0) %>%
        pull("splicejnc")
    nsj <- sum(gi$splicejnc)

    gi$idx <- NULL
    gi$intra <- NULL
    gi$gap <- NULL
    gi$same_strand <- NULL


    message("Chimeras to test against splice junctions: ", length(gi))
    message("SJ entries  : ", nsj, ": ", round(nsj / length(gi) * 100, 3), " % of all")
    message("SJ entries / single chr entries: ", nsj, "/", length(intragi), " : ", round(nsj / length(intragi) * 100, 2), " % of single chromosome chimeras")

    return(gi)
}
