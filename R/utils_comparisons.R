#' Compare multiple RNA-RNA interactions sets
#'
#' Combines all interaction into single superset by clustering & collapsing.
#' Then compares every input entry with the superset. Overlaps between superset and
#' inputs are recorded in a table as 0/1
#'
#' @param gi_samples_list anmes list with the `GInteractions` entries list('sample1'=gi1,'sample2'='gi2)
#' @param min_ratio If the overlap-to-span ratio for either arm (A or B) for pair
#' of chimeric reads is less than \code{min_arm_ratio}, then the total overlap for this pair is set to zero.
#' Relevant to comparison of superset vs individual samples
#' @param minoverlap Parameter for read clustering to create a superset. Minimum required overlap to for either arm (A or B) for pair of entries.
#' @param maxgap Parameter for read clustering. Minimum required shift between start and end coordinates of arms for pair of overlapping entries..
#' If the shift is longer than \code{max_gap} for either arm, then total read overlap between those reads is zero.
#' @param gi_superset Optional. Superset defining the space (all) of the interactions, against which inputs from the list will be compared.
#' @param niter Internal parameter for debugging. Number of cluster& collapse iterations to find superset
#' @param anno_gr Optional. `Granges` to annotate superset.
#' @return dataframe recodding the overlaps between samples and supeset
#' @export
#' @examples
#' # Create test set of RNA interactions
#' chrom <- "chr1"
#' start1 <- c(1, 11, 21, 31, 41, 51, 61, 71, 81, 91)
#' end1 <- start1 + 9
#' start2 <- c(101, 111, 121, 131, 141, 151, 161, 171, 181, 191)
#' end2 <- start2 + 9
#'
#' anchor1 <- GRanges(seqnames = chrom, ranges = IRanges(start = start1, end = end1))
#' anchor2 <- GRanges(seqnames = chrom, ranges = IRanges(start = start2, end = end2))
#'
#' interaction <- GInteractions(anchor1, anchor2)
#'
#' # Ensure some overlaps
#' n <- length(interaction)
#' group_size <- ceiling(n / 2)
#' group_indices1 <- sort(sample(seq_len(n), group_size))
#' group_indices2 <- sort(sample(seq_len(n), group_size))
#' group_indices3 <- sort(sample(seq_len(n), group_size))
#'
#' # Create separate GInteractions objects for each group
#' group1 <- interaction[group_indices1]
#' group2 <- interaction[group_indices2]
#' group3 <- interaction[group_indices3]
#'
#' # format input and call comparison
#' a <- list("sample1" = group1, "sample2" = group2, "sample3" = group3)
#' res <- compareMultipleInteractions(a)
#' # comparison result
#' head(res$dt_upset)
#' # superset
#' res$gi_all
#' # dataframe for the Upset plot
#' res$dt_upset
compareMultipleInteractions <- function(gi_samples_list,
    min_ratio = 0.3,
    minoverlap = 5,
    maxgap = 50,
    niter = 3,
    gi_superset = NULL,
    anno_gr = NULL) {
    if (!is.null(gi_superset)) {
        message("Using provided superset of interactions")
        gi_all <- gi_superset
        gi_all$n_reads <- 1

        gi_all <- refresh_gi(gi_all)
        gi_all$duplex_id <- seq_len(length(gi_all))

        statdf <- tibble(
            "read_id" = seq_len(length(gi_all)),
            "duplex_id" = seq_len(length(gi_all))
        )
        gi_all_inter <- gi_all
        gi_all_inter$id <- seq_len(length(gi_all_inter))
    } else {
        message("Gathering superset of interactions")
        gi_all <- GInteractions()
        for (i in seq_len(length(gi_samples_list))) {
            samplename <- names(gi_samples_list[i])
            gi_all <- c(gi_all, gi_samples_list[[i]])
        }
        gi_all <- refresh_gi(gi_all)
        gi_all$duplex_id <- seq_len(length(gi_all))
        gi_all$score <- 1
        statdf <- tibble(
            "read_id" = seq_len(length(gi_all)),
            "duplex_id" = seq_len(length(gi_all))
        )

        # message("DEBUG INTER:", length(gi_all))
        # gi_all_inter = clusterDuplexGroups(gi = gi_all,
        #                     maxgap =maxgap, minoverlap =  minoverlap,
        #                     min_arm_ratio = min_ratio)
        # gi_all_inter = collapse_duplex_groups(gi_all_inter,
        #                                       return_unclustered = TRUE)

        gi_all_inter <- collapseSimilarChimeras(
            gi = gi_all, niter = niter, minoverlap = minoverlap,
            maxgap = maxgap, read_stats_df = statdf
        )$gi_updated
        gi_all_inter$id <- seq_len(length(gi_all_inter))
    }


    if (!is.null(anno_gr)) {
        gi_all_inter <- annotateGI(gi_all_inter, anno_gr)
        gi_all_inter <- .annotateCisTrans(gi_all_inter)
    }


    dt_upset <- (data.frame(matrix(
        ncol = length((gi_samples_list)) + 1,
        nrow = length(gi_all_inter),
        data = 0
    )))
    colnames(dt_upset) <- c("id", names(gi_samples_list))
    dt_upset$id <- gi_all_inter$id
    if (!is.null(gi_all_inter$cis)) {
        dt_upset$cis <- gi_all_inter$cis
    }
    dt_upset$n_reads <- gi_all_inter$n_reads
    for (i in seq_len(length(gi_samples_list))) {
        samplename2 <- names(gi_samples_list[i])
        message("Superset vs ", samplename2)
        first <- gi_all_inter
        second <- assign_name_to_gi(gi_samples_list[[i]], sample_name = samplename2)
        # unset maxgap because we compare against less-defined superset
        ovldt <- computeGIPairOverlaps(first, second,
            minoverlap = minoverlap,
            maxgap = 100
        )

        # weak_overlaps_n <- nrow(ovldt %>%
        #     dplyr::filter(!(ratio.A >= min_ratio & ratio.B >= min_ratio)))

        # message("Weak overlaps : ", weak_overlaps_n)
        # message("n_overlaps b : ",nrow(ovldt) )
        # ovldt <- ovldt %>%
        #     dplyr::filter(ratio.A >= min_ratio & ratio.B >= min_ratio)

        # message("n_overlaps f : ",nrow(ovldt) )
        if (nrow(ovldt) == 0) {
            message("no overlap")
            next
        }
        ovl_interaction_ids <- unique(ovldt$id.1)
        dt_upset[dt_upset$id %in% ovl_interaction_ids, ][samplename2] <- 1
    }
    res <- list()
    res$gi_all <- gi_all_inter
    res$dt_upset <- dt_upset
    return(res)
}

computeGIPairOverlaps <- function(gi1, gi2, minoverlap = 10,
    maxgap = 10) {
    # this is to unify levels
    # otherwise A-B B-A matches would be messed between query an dsubject

    mcols(gi1) <- NULL
    mcols(gi2) <- NULL
    runi <- unify_gi_levels(gi1, gi2)
    gi1 <- runi$gi1
    gi2 <- runi$gi2
    fo <- findOverlaps(gi1, gi2,
        ignore.strand = FALSE,
        type = "equal", use.region = "both", maxgap = maxgap,
        minoverlap = minoverlap, select = "all"
    )

    dt <- tibble(
        id.1 = queryHits(fo),
        id.2 = subjectHits(fo),
        A_span = width(punion(gi1@regions[gi1[queryHits(fo)]@anchor1],
            gi2@regions[gi2[subjectHits(fo)]@anchor1],
            fill.gap = TRUE
        )),
        B_span = width(punion(gi1@regions[gi1[queryHits(fo)]@anchor2],
            gi2@regions[gi2[subjectHits(fo)]@anchor2],
            fill.gap = TRUE
        )),
        A_ovl = width(pintersect(
            gi1@regions[gi1[queryHits(fo)]@anchor1],
            gi2@regions[gi2[subjectHits(fo)]@anchor1]
        )),
        B_ovl = width(pintersect(
            gi1@regions[gi1[queryHits(fo)]@anchor2],
            gi2@regions[gi2[subjectHits(fo)]@anchor2]
        ))
    ) %>% mutate(ratio.A = (A_ovl / A_span), ratio.B = (B_ovl) / (B_span))
    gc()
    return(dt)
}
