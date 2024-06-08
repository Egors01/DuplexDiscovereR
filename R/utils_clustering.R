#' Helper function to add ids to the duplex groups missed during global clustering
#'
#' @description
#' Check if there are a temporary duplex records with \code{duplex_id}, which consist of more than one read  \code{n_reads > 1} , but
#' does not have assigned any \code{dg_id} as the duplex group (DG) index. Creates new  \code{dg_id} if \code{n_reads > 1}
#' @details
#' Meant to be used in the situations when previous collapsing steps merged two or more reads to the temporary DG with \code{duplex_id}, but
#' global clustering has not identified any overlap between this temporary group and other duplexes, resulting in undefined \code{dg_id}.
#' This function looks up for these cases and creates new \code{dg_id} for temporary DGs, marking them as the final DGs.
#' New \code{dg_id} values are unique and allocated sequentially after the  maximum value of \code{dg_id}
#' @param gi_input \pkg{GInteractions} with the  \code{dg_id},  \code{duplex_id} and  \code{n_reads} column
#' @return
#' \pkg{GInteractions} object with new \code{dg_id} for rows with  \code{n_reads > 1}
add_dg_ids_for_noempty_duplexes <- function(gi_input) {
    # add new cluster ids
    col_check <- c("duplex_id", "dg_id", "n_reads")
    if (!all(col_check %in% colnames(mcols(gi_input)))) {
        ermsg <- paste0(col_check, sep = " ")
        stop("One of the ", ermsg, " columns is not in the input gi ")
    }

    n_add <- length(gi_input[is.na(gi_input$dg_id) & gi_input$n_reads > 1])
    if (n_add == 0) {
        message("No collapsed duplexes with n_reads>1 and without dg_id found")
        return(gi_input)
    } else {
        # max_dg_id = max(gi_input$dg_id,na_rm=TRUE) # returns NA for some reason
        max_dg_id <- max(gi_input$dg_id[!is.na(gi_input$dg_id)])

        gi_keep <- gi_input[!(is.na(gi_input$dg_id) & gi_input$n_reads > 1)]

        gi_change <- gi_input[is.na(gi_input$dg_id) & gi_input$n_reads > 1]
        df_add <- makeDfFromGi(gi_change) %>%
            mutate(dg_id = max_dg_id + row_number())
        mcols(gi_change) <- df_add

        # message(" Changing ",n_add, "records from NA dg_id to",paste0(df_add$dg_id,sep=' '))
        gi_res <- c(gi_keep, gi_change)
        return(gi_res)
    }
}


#' Call clustering multiple times to collapse similar reads into duplex groups
#'
#' Function calls clustering algorithm several times and collapses highly similar
#' reads to the temporary duplex groups (DGs).
#'
#' @details
#' Calling this procedure before global read clustering
#' substantially reduces time required for calling DGs.
#' Collapsed duplex groups are aggregated only from the reads which are shifted
#'  by only a few nucleotides from each other.
#' These DGs are temporary until full library clustering is called.
#' To keep track of the mapping of the temprary DGs to the input, dedicated
#' dataframe is returned. The 'duplex_id' column will be added or updated as
#' identifier for the temporary duplex group.
#' The number of reads under single 'duplex_id' is recorded in the 'n_reads' fields
#' @param gi  `GInteractions` object
#' @param read_stats_df `tibble` with the mapping 'read_id' and 'duplex_id' fields
#'  'read_id' refers to the unique read, 'duplex_id' refers to the entry collapsed
#'  identical reads i.e two identical reads will will correspond to two unique read_id and
#'  the single duplex_id with n_reads=2
#' @param mgap  Maximum relative shift between the overlapping read arms
#' @param niter Number of times clustering will be called
#' @param minovl Minimum required overlap between either read arm
#'
#' @return a list with the following keys
#' \describe{
#'   \item{gi_updated}{ `GInteractions` object with both collapsed duplex groups
#'   and not-collapsed unchanged reads}
#'   \item{stats_df}{ `tibble` With the mapping from the unique read -
#'    with the the infromation about time and memory reaquired for the function
#'    call}
#' }
collapse_similar_reads <- function(
        gi, read_stats_df,
        mgap = 5,
        niter = 2,
        minovl = 10) {
    message("--- Collapsing of the reads shifted by <= ", mgap, " nt ---")
    gi_base <- gi
    gi_base$dg_id <- NULL
    stats_df <- read_stats_df
    read_stats_df_upd <- read_stats_df
    for (i in seq_len(niter)) {
        ovl_df <- compute_gi_self_overlaps(gi_base, maxgap = mgap, minovl = minovl)
        message("----iter ", i, "-----")
        message(
            "Connectivity of the duplex graph: ", get_connectivity(gi_base, ovl_df),
            " n_edges = ", nrow(ovl_df), " n_nodes = ", length(gi_base)
        )
        if (is_empty(ovl_df) | nrow(ovl_df) < 3) {
            message("nothing to collapse on iter", i)
            break
        }

        gi_clustered <- clusterDuplexGroups(gi = gi_base, graphdf = ovl_df, maxgap = mgap)

        grc <- collapse_duplex_groups(gi_clustered,
            return_unclustered = TRUE,
            return_collapsed = FALSE
        )
        message("N nodes after collapse: ", length(grc))

        # update duplex_ids
        orig_mcols <- as_tibble(data.frame(mcols(gi_clustered)[c("duplex_id", "dg_id")]))
        new_mcols <- as_tibble(data.frame(mcols(grc)[c("duplex_id", "dg_id")]))
        fill_mcols <- orig_mcols %>%
            dplyr::filter(!is.na(dg_id)) %>%
            group_by(dg_id) %>%
            dplyr::summarize(
                new_duplex_id = dplyr::first(duplex_id)
            )


        # update in the gi object
        new_collapsed_stats <- left_join(new_mcols, fill_mcols, by = "dg_id") %>%
            mutate(
                new_duplex_id = ifelse(!is.na(duplex_id), duplex_id, new_duplex_id)
            )
        grc$duplex_id <- new_collapsed_stats %>% pull(new_duplex_id)
        if (any(is.na(grc$duplex_id))) {
            stop("After collapsing: Index column 'duplex_id' contains undefined variables")
        }
        # update the stats dataframe parallel to the original gi object
        new_original_stats <- left_join(orig_mcols, fill_mcols, by = c("dg_id")) %>%
            mutate(
                old_duplex_id = duplex_id,
                new_duplex_id = ifelse(is.na(dg_id), old_duplex_id, new_duplex_id)
            ) %>%
            dplyr::select(c(new_duplex_id, old_duplex_id, dg_id))

        #
        # read_stats_df_upd =  bind_cols(
        #   stats_df %>%  dplyr::select(-duplex_id),
        #   new_original_stats %>%  dplyr::select(-c(old_duplex_id,dg_id))) %>%
        #   dplyr::rename(duplex_id=new_duplex_id)

        read_stats_df_upd <- left_join(stats_df, new_original_stats,
            by = c("duplex_id" = "old_duplex_id")
        ) %>%
            mutate(duplex_id = if_else(is.na(dg_id), duplex_id, new_duplex_id)) %>%
            dplyr::select(-c(new_duplex_id, dg_id)) %>%
            dplyr::distinct(read_id, .keep_all = TRUE)

        # update for next round or exit
        grc$dg_id <- NULL
        gi_base <- grc
        stats_df <- read_stats_df_upd
    }
    res <- list()
    res$gi_updated <- gi_base
    res$stats_df <- read_stats_df_upd

    return(res)
}


#' Accessor for mapping between temporary and final cluster ids
#' @keywords internal
#' @param gi
#'
#' @return tibble
dg_id_to_duplex_id <- function(gi) {
    future_clusters <- as_tibble(data.frame(gi))
    return(future_clusters %>% dplyr::select(duplex_id, dg_id))
}

#' Collapses identical interactions
#' @description
#' Two entries (reads) are considered identical if they share start, end, strand and score vales
#' Identical entries are collapsed into the single one.
#' @details
#'
#' Adds columns to the collapsed object
#' duplex_id (int) unique record id
#' n_reads (int) number of entries collapsed
#'
#' @param gi GInteractions(mode='strict') object with chromA, strandA, startA,
#' endA, chromB, strandB, startB, endB, score columns
#' Optionally cigar_alnA, cigar_alnB columns are also considered for collapsing
#' 'read_id' column used as the index in the initial objects. Created, if not exists
#' @return result_list object with keys
#'' gi_collapsed':   New collapsed GInteraction object
#'' stats_df': tibble with the mapping of the original entries to the new duplex_id
#' @export
#' @examples
#' # load data
#' data("RNADuplexesSmallGI")
#' res_collapse <- collapse_identical_reads(SampleSmallGI)
#' gi_new <- res_collapse[["gi_collapsed"]]
#' # keeps the mapping of the colapsed object to new
#' read_stats_df <- res_collapse[["stats_df"]]
collapse_identical_reads <- function(gi) {
    if (is(gi, "StrictGInteractions")) {
        out_gi <- TRUE
        gi_dt <- makeDfFromGi(gi)
    } else {
        out_gi <- FALSE
        message("Input type is not StrictGInteractions, will try to work with the dataframe")
        message("Expected columns: chromA, strandA, startA, endA, chromB, strandB, startB, endB, score")
        message("Optional: cigar_alnA, cigar_alnB")
        gi_dt <- gi
    }

    if (!("read_id" %in% colnames(gi_dt))) {
        gi_dt$read_id <- seq_len(nrow(gi_dt))
        message("Created 'read_id' column as the unique index in  the provided object")
    }

    if ("cigar_alnA" %in% colnames(gi)) {
        dt1 <- gi_dt %>%
            tidyr::unite("grp", c(
                chromA, strandA, startA, endA,
                chromB, strandB, startB, endB, cigar_alnA, cigar_alnB, score
            ), remove = FALSE) %>%
            mutate(num = digest2int(grp))
    } else {
        dt1 <- gi_dt %>%
            tidyr::unite("grp", c(
                chromA, strandA, startA, endA,
                chromB, strandB, startB, endB, score
            ), remove = FALSE) %>%
            mutate(num = digest2int(grp))
    }

    dt1$n_reads <- NULL

    readcts <- dt1 %>%
        dplyr::filter(duplicated(num)) %>%
        dplyr::select(num)
    readcts <- table(readcts$num)
    readcts <- tibble("num" = as.double(names(readcts)), "n_reads" = as.vector(readcts) + 1)

    # dt_or = dt1

    # dt1 = dt1 %>% dplyr::distinct(num,.keep_all = TRUE)
    dt1 <- left_join(dt1, readcts, by = "num") %>%
        mutate(n_reads = tidyr::replace_na(n_reads, 1))

    read2duplex_map <- dt1 %>% dplyr::select(read_id, num)

    dt1 <- dt1 %>% dplyr::distinct(num, .keep_all = TRUE)
    dt1 <- dt1 %>%
        mutate(duplex_id = seq_len(nrow(dt1)))

    read2duplex_map <- left_join(read2duplex_map, dt1[, c("num", "duplex_id", "n_reads")], by = "num") %>%
        dplyr::rename(n_reads_collapsed = n_reads) %>%
        dplyr::select(-c(num))
    # dt1 =dt1 %>%
    #   dplyr::select(-c(grp,num,read_id))
    #
    dt1 <- dt1 %>%
        dplyr::select(c(
            chromA, strandA, startA, endA,
            chromB, strandB, startB, endB, duplex_id, n_reads, score
        ))


    dt1 <- dt1 %>%
        relocate(duplex_id, .before = score)
    message("Duplicated   :  ", round((1 - nrow(dt1) / nrow(gi_dt)) * 100, 2), "% of initial")
    message("Initial size :  ", nrow(gi_dt))
    message("New size     :  ", nrow(dt1))

    if (out_gi == TRUE) {
        dt1 <- makeGiFromDf(dt1)
    }

    res <- list()
    res$gi_collapsed <- dt1
    res$stats_df <- read2duplex_map
    return(res)
}




#' Collapses identical interactions
#' @description
#' Two entries (reads) are considered identical if they share start, end, strand and score vales
#' Identical entries are collapsed into the single one.
#' @details
#' Adds columns to the collapsed object
#' duplex_id (int) unique record id
#' n_reads (int) number of entries collapsed
#' @param gi GInteractions(mode='strict') object with chromA, strandA, startA, endA, chromB, strandB, startB, endB, score columns
#' Optionally cigar_alnA, cigar_alnB columns are also considered for collapsing
#' 'read_id' column used as the index in the initial objects. Created, if not exists
#' @return result_list object with keys
#'' gi_collapsed':   New collapsed GInteraction object
#'' stats_df': tibble with the mapping of the original entries to the new duplex_id
#' @export
#' @examples
#' # load data
#' data("RNADuplexesSmallGI")
#' res_collapse <- collapse_identical_reads(SampleSmallGI)
#' gi_new <- res_collapse[["gi_collapsed"]]
#' read_stats_df <- res_collapse[["stats_df"]]
collapse_identical_reads <- function(gi) {
    if (is(gi, "StrictGInteractions")) {
        out_gi <- TRUE
        gi_dt <- makeDfFromGi(gi)
    } else {
        out_gi <- FALSE
        message("Input type is not StrictGInteractions, will try to work with the dataframe")
        message("Expected columns: chromA, strandA, startA, endA, chromB, strandB, startB, endB, score")
        message("Optional: cigar_alnA, cigar_alnB")
        gi_dt <- gi
    }

    if (!("read_id" %in% colnames(gi_dt))) {
        gi_dt$read_id <- seq_len(nrow(gi_dt))
        message("Created 'read_id' column as the unique index in  the provided object")
    }

    if ("cigar_alnA" %in% colnames(gi)) {
        dt1 <- gi_dt %>%
            tidyr::unite("grp", c(
                chromA, strandA, startA, endA,
                chromB, strandB, startB, endB, cigar_alnA, cigar_alnB, score
            ), remove = FALSE) %>%
            mutate(num = digest2int(grp))
    } else {
        dt1 <- gi_dt %>%
            tidyr::unite("grp", c(
                chromA, strandA, startA, endA,
                chromB, strandB, startB, endB, score
            ), remove = FALSE) %>%
            mutate(num = digest2int(grp))
    }

    dt1$n_reads <- NULL

    readcts <- dt1 %>%
        dplyr::filter(duplicated(num)) %>%
        dplyr::select(num)
    readcts <- table(readcts$num)
    readcts <- tibble("num" = as.double(names(readcts)), "n_reads" = as.vector(readcts) + 1)

    # dt_or = dt1

    # dt1 = dt1 %>% dplyr::distinct(num,.keep_all = TRUE)
    dt1 <- left_join(dt1, readcts, by = "num") %>%
        mutate(n_reads = tidyr::replace_na(n_reads, 1))

    read2duplex_map <- dt1 %>% dplyr::select(read_id, num)

    dt1 <- dt1 %>% dplyr::distinct(num, .keep_all = TRUE)
    dt1 <- dt1 %>%
        mutate(duplex_id = seq_len(nrow(dt1)))

    read2duplex_map <- left_join(read2duplex_map, dt1[, c("num", "duplex_id", "n_reads")], by = "num") %>%
        dplyr::rename(n_reads_collapsed = n_reads) %>%
        dplyr::select(-c(num))
    # dt1 =dt1 %>%
    #   dplyr::select(-c(grp,num,read_id))
    #
    dt1 <- dt1 %>%
        dplyr::select(c(
            chromA, strandA, startA, endA,
            chromB, strandB, startB, endB, duplex_id, n_reads, score
        ))


    dt1 <- dt1 %>%
        relocate(duplex_id, .before = score)
    message("Duplicated   :  ", round((1 - nrow(dt1) / nrow(gi_dt)) * 100, 2),
            "% of initial")
    message("Initial size :  ", nrow(gi_dt))
    message("New size     :  ", nrow(dt1))

    if (out_gi == TRUE) {
        dt1 <- makeGiFromDf(dt1)
    }

    res <- list()
    res$gi_collapsed <- dt1
    res$stats_df <- read2duplex_map
    return(res)
}

get_connectivity <- function(gi, dt_conn) {
    edg <- nrow(dt_conn)
    vert <- length(gi)
    connectivity <- round(edg / vert, 4)
    return(connectivity)
}

#' Find overlaps between entries in `GInteractions`
#'
#' Utility function to find overlapping reads in the input. Removes
#' self-hits. Computes overlap/span ratios for
#' each interaction arm. Sum of the scores is recorded in 'weight' field
#'
#' @param tr input gi object
#' @param id_column column which use for using as ids for entries
#' @param maxgap parameter for call of [InteractionSet::findOverlaps()]
#' @param minovl  parameter for call [InteractionSet::findOverlaps()]
#'
#' @return dataframe with indexes of pairwise overlapsin input and
#' columns for span, overlap, ratios of either read arm
#' @export
#' @examples
#' data("RNADuplexesSmallGI")
#' compute_gi_self_overlaps(SampleSmallGI)
compute_gi_self_overlaps <- function(tr, id_column = "duplex_id", maxgap = 40, minovl = 10) {
    if (id_column %in% colnames(mcols(tr))) {
        # message("Using ",id_column," as id for computing overlaps")
        ids <- mcols(tr)[, id_column]
        mcols(tr) <- NULL
        tr$idcol <- ids
    } else {
        # message("Using index as id for computing overlaps")
        mcols(tr) <- NULL
        tr$idcol <- seq_len(length(tr))
        id_column <- "index"
    }
    id_names <- str_c(id_column, c(1, 2), sep = ".")

    fo <- findOverlaps(tr, ignore.strand = FALSE, type = "equal", maxgap = maxgap, minoverlap = minovl)
    fo <- fo[queryHits(fo) < subjectHits(fo)]
    if (length(fo) == 0) {
        dt <- tibble()
        message("No overlap found")
        return(dt)
    }
    dt <- tibble(
        idcol.1 = tr[queryHits(fo)]$idcol,
        idcol.2 = tr[subjectHits(fo)]$idcol,
        A_span = width(punion(
            tr@regions[tr[queryHits(fo)]@anchor1],
            tr@regions[tr[subjectHits(fo)]@anchor1]
        )),
        B_span = width(punion(
            tr@regions[tr[queryHits(fo)]@anchor2],
            tr@regions[tr[subjectHits(fo)]@anchor2]
        )),
        A_ovl = width(pintersect(
            tr@regions[tr[queryHits(fo)]@anchor1],
            tr@regions[tr[subjectHits(fo)]@anchor1]
        )),
        B_ovl = width(pintersect(
            tr@regions[tr[queryHits(fo)]@anchor2],
            tr@regions[tr[subjectHits(fo)]@anchor2]
        ))
    ) %>%
        mutate(
            ratio.A = (A_ovl / A_span),
            ratio.B = (B_ovl) / (B_span),
            weight = ratio.A + ratio.B
        ) %>%
        dplyr::rename(
            !!quo_name(id_names[1]) := idcol.1,
            !!quo_name(id_names[2]) := idcol.2
        )

    return(dt)
}
