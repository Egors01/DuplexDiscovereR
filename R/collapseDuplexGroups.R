#' Collapse the reads into the duplex groups after clustering
#'
#' @description
#' Collapse each interaction in the input to the duplex group based on the pre-computed dg_id
#' @param gi `GInteractions` with the 'dg_id' metadata field
#' @param return_unclustered add unclustered reads to output
#' @param return_collapsed add duplex groups, which were created as temporary
#' with n_reads > 1 but was not clustered to the DG golabally. This parameter
#' is used internally and should be kept default in most situations.
#' @param keep_meta whether to keep metadata, which only unclustered reads have, in case of a mixed output
#'
#' @details
#' 'dg_id' is used as the identifier for the duplex group
#' Reads belonging to the same duplex group are collapsed into a single entry
#' with start and end are set as min() and max() coordinate of the reads in within the duplex group.
#' The 'score' column is averaged across the duplex group reads is calculated and
#' put as the 'score' for the collapsed duplex group
#' Behavior in case 'dg_id'  = NA: Option '\code{return_unclustered}' -
#'  whether unclustered reads with should be added to the output gi
#' \describe{
#'   \item{return_unclustered == FALSE}{Interaction is not returned in the output. Default. }
#'   \item{return_unclustered == TRUE}{Interaction is returned in the output,
#'    output is mixed duplex groups and individual reads}
#'   }
#'  Internally used argument
#'  #' \describe{
#'   \item{return_collapsed == FALSE}{In case interaction already collapsed and n_read > 1,
#'   interaction will not be returned as duplex group }
#'   \item{return_collapsed == TRUE}{In case interaction has n_read > 1, interaction will be treated
#'    as duplex group }
#'   }
#' @importFrom dplyr first
#' @export
#' @return
#' `GInteractions` object with collapsed duplex groups
#' @examples
#' # load example of clustered data
#' data("RNADuplexesSampleData")
#' # some reads assigned to DG, some are not
#' table(is.na(RNADuplexSampleGI$dg_id))
#' # Return only DGs
#' gicollapsed <- collapse_duplex_groups(RNADuplexSampleGI, return_unclustered = FALSE)
#' # Return DGs and unclustered reads as well
#' gimixed <- collapse_duplex_groups(RNADuplexSampleGI, return_unclustered = TRUE)
#'
#' # load small sample GInteractions and process it manually
#' data("RNADuplexesSmallGI")
#' # First, collapse duplicated reads. This adds n_reads and duplex ids
#' ginodup <- collapseIdenticalReads(SampleSmallGI)$gi_collapsed
#' # Second, run clustering, get DG ids
#' ginodup <- clusterDuplexGroups(ginodup)
#' # Return all DGs result in n=3 DGS, one of them formed by
#' # identical duplicated alignments
#' collapse_duplex_groups(ginodup, return_collapsed = TRUE)
#' # Return DGs, but drop duplicated returns n=2 DGs
#' collapse_duplex_groups(ginodup, return_collapsed = FALSE)
collapse_duplex_groups <- function(
        gi, return_unclustered = FALSE,
        return_collapsed = TRUE,
        keep_meta = TRUE) {
    # check n_reads field
    if (!is.null(gi$n_reads)) {
        gi$n_reads <- as_tibble(mcols(gi)["n_reads"]) %>%
            mutate(n_reads = if_else(!is.na(n_reads), n_reads, 1)) %>%
            pull(n_reads)
    } else {
        gi$n_reads <- 1
    }

    if (return_collapsed) {
        to_add_n <- length(gi[is.na(gi$dg_id) & gi$n_reads > 1])
        if (to_add_n != 0) {
            message("Found temporary duplex groups without assigned dg_id", to_add_n)
            gi <- .addDGidsForTmpDGs(gi)

            # Consider calling  .addDGidsForTmpDGs(gi) to add them into final DGs")
        }
    }

    future_clusters <- as_tibble(data.frame(gi))

    future_clusters <- future_clusters %>%
        dplyr::filter(!is.na(dg_id)) %>%
        group_by(dg_id) %>%
        summarise(
            chromA = dplyr::first(seqnames1),
            startA = min(start1),
            endA = max(end1),
            strandA = dplyr::first(strand1),
            chromB = dplyr::first(seqnames2),
            startB = min(start2),
            endB = max(end2),
            strandB = dplyr::first(strand2),
            n_reads = sum(n_reads),
            score = mean(score)
        ) %>%
        relocate(dg_id, .after = n_reads)

    res_gi <- makeGiFromDf(future_clusters)
    if (return_unclustered) {
        unclustered <- gi[is.na(gi$dg_id)]
        res_gi <- c(res_gi, unclustered)
    }

    return(res_gi)
}
