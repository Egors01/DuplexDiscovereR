#' Get colnames for expected data types
#'
#' @param nameset name of the table
#' @keywords internal
#' @return character vector
get_colnames_and_types_for_input <- function(nameset) {
    switch(nameset,
        "STAR_COLNAMES19" = c(
            "chr_donorA", "brkpt_donorA", "strand_donorA", "chr_acceptorB",
            "brkpt_acceptorB", "strand_acceptorB", "junction_type", "repeat_left_lenA",
            "repeat_right_lenB", "read_name", "start_alnA", "cigar_alnA",
            "start_alnB", "cigar_alnB", "num_chim_aln", "max_poss_aln_score",
            "non_chim_aln_score", "this_chim_aln_score", "bestall_chim_aln_score"
        ),
        "STAR_COLNAMES14" = c(
            "chr_donorA", "brkpt_donorA", "strand_donorA", "chr_acceptorB",
            "brkpt_acceptorB", "strand_acceptorB", "junction_type", "repeat_left_lenA",
            "repeat_right_lenB", "read_name", "start_alnA", "cigar_alnA",
            "start_alnB", "cigar_alnB"
        ),
        "STAR_COLTYPES19" = c(
            "character", "numeric", "character", "character", "numeric",
            "character", "numeric", "numeric", "numeric", "character",
            "numeric", "character", "numeric", "character", "numeric",
            "numeric", "numeric", "numeric", "numeric"
        ),
        "STAR_COLTYPES14" = c(
            "character", "numeric", "character", "character", "numeric",
            "character", "numeric", "numeric", "numeric", "character",
            "numeric", "character", "numeric", "character"
        ),
        "BEDPE_COLNAMES" = c(
            "chromA", "startA", "endA", "chromB",
            "startB", "endB", "readname", "flag", "strandA", "strandB"
        ),
        "BEDPE_COLTYPES" = c(
            "character", "numeric", "numeric", "character",
            "numeric", "numeric", "character", "numeric", "character", "character"
        ),
        stop("Invalid variable name")
    )
}


#' Count the length of the key type in CIGAR string
#' @description
#'
#' Takes CIGAR operands i.e M,N,S and sums the associated blocks length
#' It is vectorized. i.e supports vector with CIGAR strings
#' @param strings CIGAR string vector
#' @param s CIGAR operands
#'
#' @return vector with length values
#' @export
#'
#' @examples
#' # From a vector
#' get_char_count_cigar(c("4S18M22S", "25S26M"), "S")
#' get_char_count_cigar(c("18M22S", "20M20S"), "M")
get_char_count_cigar <- function(strings, s) {
    reg_tmp <- "\\d+(?=V)"
    reg <- str_replace(reg_tmp, "V", s)
    charcount <- rowSums(as_tibble(data.frame(
        str_extract_all(strings, reg,
            simplify = TRUE
        ),
        fix.empty.names = TRUE
    )) %>%
        mutate_all((as.numeric)) %>%
        replace(is.na(.), 0) %>% as.matrix())

    return(charcount)
}



#' Convert Dataframe to GInteractions
#'
#' @description
#' Converts dataframe-like object to the `GInteractions`.
#' @details
#' arms will be consistent between different objects of same reference
#' Following columns are looked up in input dataframe  to parse region coordinates:
#' c("chromA','startA','endA','strandA',"chromB",'startB','endB','strandB')
#' `GInteractions(mode='strict')` is enforced, to ensure that the order of the regions
#' Extra columns are stored as metadata fields
#' @param df dataframe-like object. Should be convertable to tibble::tibble()
#' @return GInteractions(mode='strict')
#' @seealso [DuplexDiscovereR::makeDfFromGi()]
#' @export
#' @examples
#' # load example GInteractions
#' data(RNADuplexesSmallGI)
#'
#' converted_to_df <- makeDfFromGi(SampleSmallGI)
#' converted_to_gi <- makeGiFromDf(converted_to_df)
makeGiFromDf <- function(df) {
    dtt <- tibble(df)
    # requires "chromA","startA","endA","strandA",
    # "chromB","startB","endB","strandB" columns
    df_a <- dtt %>%
        dplyr::select(chromA, startA, endA, strandA) %>%
        purrr::set_names(c("chrom", "start", "end", "strand"))

    df_b <- dtt %>%
        dplyr::select(chromB, startB, endB, strandB) %>%
        purrr::set_names(c("chrom", "start", "end", "strand"))
    meta <- dtt %>%
        dplyr::select(!dplyr::one_of(c(
            "chromA", "startA", "endA", "strandA",
            "chromB", "startB", "endB", "strandB"
        )))
    a_gr <- GRanges(df_a)
    b_gr <- GRanges(df_b)
    gi <- GInteractions(a_gr, b_gr, mode = "strict")
    mcols(gi) <- data.frame(meta)
    return(gi)
}

#' Convert GInteractions to tibble
#'
#' @description
#' Converts `GInteractions` to tibble, preserves metadata
#' @details
#' Following naming conventions is used  for region coordinates:
#' `c('chromA','startA','endA','strandA',
#' 'chromB','startB','endB','strandB')`
#' @param gi GInteracttions
#' @return tibble preserving metadata columns
#' @seealso [DuplexDiscovereR::makeGiFromDf()]
#' @export
#' @examples
#' data(RNADuplexesSmallGI)
#' converted_to_df <- makeDfFromGi(SampleSmallGI)
#' converted_to_gi <- makeGiFromDf(converted_to_df)
makeDfFromGi <- function(gi) {
    meta <- as_tibble(mcols(gi))
    mcols(gi) <- NULL
    gidf <- as_tibble(data.frame(gi)) %>%
        dplyr::rename(
            chromA = seqnames1,
            startA = start1,
            endA = end1,
            strandA = strand1,
            chromB = seqnames2,
            startB = start2,
            endB = end2,
            strandB = strand2,
        ) %>%
        dplyr::select(-c(width1, width2))
    gidf <- bind_cols(gidf, meta)
    return(gidf)
}

#' Get chimeric junctions
#'
#' Returns chimeric junction defined as range distance between the end and the
#' start of the first and second range respectively
#'
#' @details
#' If the pair of interacting ranges is not on the same strand and chromosome,
#' returns error
#' @keywords internal
#' @returns Granges object with the
get_chimeric_junctions_onestrand <- function(gi_intra) {
    gi_intra$intra <- pairdist(gi_intra, type = "intra")
    if (!all(gi_intra$intra == TRUE)) {
        stop("Cannot extract chimeric junction. Not all ranges are on the same chromosome/strand")
    }
    gr_chim <- GRanges(data.frame(
        "chrom" = seqnames(gi_intra@regions[gi_intra@anchor2]),
        "start" = end(gi_intra@regions[gi_intra@anchor1]),
        "end" = start(gi_intra@regions[gi_intra@anchor2]),
        "strand" = strand(gi_intra@regions[gi_intra@anchor2])
    ))
    return(gr_chim)
}


#' Refresh the `GInteractions` object
#'
#' Sub-setting the `GInteractions` object does not reduce its `ranges` container
#' For some applications, to save memory, we can safely reduce the size of the
#' object by re-creating it.
#' Also, it can be used to ensure the 'strict' mode of the regions in `ranges`
#'
#' @param gi `GInteractions` object
#' @keywords internal
#' @return `GInteractions` object with new ranges attribute
refresh_gi <- function(gi) {
    a1 <- gi@anchor1
    a2 <- gi@anchor2
    gi_collapsed <- GInteractions(gi@regions[a1], gi@regions[a2], mode = "strict")
    mcols(gi_collapsed) <- mcols(gi)
    return(gi_collapsed)
}
#' Subset the `GInteractions` object to single interaction
#'
#' Sub-setting the `GInteractions` object does not reduce its `ranges` container
#' This function selects `Ginteraction` by index and reduces the `ranges`
#' @param gi `GInteractions` object
#' @keywords internal
#' @return `GInteractions` with the range atribure reduced to single interaction

subset_gi <- function(gi, k) {
    gi <- gi_trans$SPLASH_RA_1
    a1 <- gi[k]@anchor1
    a2 <- gi[k]@anchor2
    gi_short <- GInteractions(1, 2, gi@regions[c(a1, a2)], mode = "strict")
    mcols(gi_short) <- mcols(gi[k])
    return(gi_short)
}

#' Convert `GInteractions` object to `Granges`
#'
#' Creates the 'long' `GRanges` by stacking the A and B arms one 'on top' of the other.
#' Adds `id` and  `group` fields as indicators of original index and interaction
#' arm (A- left arm, B- right arm)
#' @param gi `GInteractions`
#' @export
#' @return `GRanges` twice the length of the input
#' @examples
#' data("RNADuplexesSmallGI")
#' convert_gi_to_ranges(SampleSmallGI)
convert_gi_to_ranges <- function(gi) {
    boxes <- c(anchors(gi)$first, anchors(gi)$second)
    groupp <- rep(seq_len(length(gi)), times = 2)
    ids <- str_c(groupp, ".", rep(c("A", "B"), each = length(gi)))
    mcols(boxes) <- rep(mcols(gi), 2)
    boxes$group <- groupp
    boxes$id <- ids
    return(boxes)
}

#' Get left arm of GInteraction
#'
#' @param gi GInteractions object
#' @keywords internal
#' @return Granges object with the left (A) region
get_arm_a <- function(gi) {
    a1 <- gi@anchor1
    a2 <- gi@anchor2
    gi_collapsed <- GInteractions(gi@regions[a1], gi@regions[a2], mode = "strict")
    return(anchors(gi_collapsed)$first)
}

#' Get right arm of GInteraction

#' @param gi GInteractions object
#' @keywords internal
#' @return Granges object with the right (B) region
get_arm_b <- function(gi) {
    a1 <- gi@anchor1
    a2 <- gi@anchor2
    gi_collapsed <- GInteractions(gi@regions[a1], gi@regions[a2], mode = "strict")
    return(anchors(gi_collapsed)$second)
}


assign_name_to_gi <- function(gi, sample_name) {
    df <- mcols(gi) %>%
        as_tibble() %>%
        mutate("sample_name" = sample_name)
    mcols(gi) <- df
    return(gi)
}

unify_gi_levels <- function(gi1, gi2) {
    gi1 <- refresh_gi(gi1)
    gi2 <- refresh_gi(gi2)
    gi1$tmpid <- 1
    gi2$tmpid <- 2
    gi3 <- c(gi1, gi2)
    gi3 <- refresh_gi(gi3)
    gi1 <- gi3[gi3$tmpid == 1]
    gi2 <- gi3[gi3$tmpid == 2]
    return(list("gi1" = gi1, "gi2" = gi2))
}
