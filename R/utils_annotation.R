#' Annotate RNA duplexes with features
#'
#' Overlays RNA duplexes with `GRanges` annotation object.
#'
#' @param gi `GInteraction` object to annotate
#' @param anno_gr `GRanges` object with the `keys` columns in the metadata
#' @param keys names of the features to use for annotation.
#' @param save_ambig in case RNA duplex overlaps multiple features of the first key,
#' mark the existense of ambiguous annotation in the fields
#' `ambig.A` and `ambig.B`. Fields `ambig_list.A` and `ambig_list.B` will be
#' store the list of overlapping features
#' Only the first filed from `keys` is checked for possible annotation
#' ambiguities.
#' @return `GInteractions` object with new fields
#' @export
#' @details
#' For each annotation feature in `keys`, i.e if keys=c(keyname1),
#' then `<keyname1>.A`, `<keyname1>.B` annotation fields will be created, containing the
#' names of overlapping features
#' If no overlap is found for the feature, then filed will have NA
#' @examples
#' data("RNADuplexesSampleData")
#' annotateGI(gi = RNADuplexSampleDGs, anno_gr = SampleGeneAnnoGR)
annotateGI <- function(
        gi, anno_gr,
        keys = c("gene_name", "gene_type", "gene_id"),
        save_ambig = TRUE) {
    orininal_columns <- colnames(mcols(gi))

    anno_gr <- anno_gr[, keys]
    anno_gr$feature_id <- seq_len(length(anno_gr))
    gi@regions$region_id <- seq_len(length(gi@regions))

    po <- findOverlapPairs(gi@regions, anno_gr)
    df <- tibble(
        "region_id" = po@first$region_id,
        "feature_id" = po@second$feature_id,
        "overlap_w" = width(pintersect(po))
    ) %>%
        group_by(region_id) %>%
        mutate(ambig = as.integer(n() > 1))
    # save ambig mapping regions (use 1st key to save ids)
    if (save_ambig) {
        ambig_save <- df %>% dplyr::filter(ambig == 1)
        ambig_save <- left_join(ambig_save,
            as_tibble(mcols(anno_gr)[c(keys[1], "feature_id")]),
            by = "feature_id"
        )
        ambig_save <- ambig_save %>%
            group_by(region_id) %>%
            dplyr::summarise(ambig_list = paste(unique(get(keys[1])), collapse = ",")) %>%
            ungroup()
    }

    df <- df %>%
        arrange(-overlap_w, .by_group = TRUE) %>%
        ungroup() %>%
        distinct(region_id, .keep_all = TRUE)
    df <- left_join(df, as_tibble(mcols(anno_gr)), by = "feature_id")

    dt <- left_join(as_tibble(mcols(gi@regions)), df, by = "region_id")

    if (save_ambig) {
        dt <- left_join(dt, ambig_save, by = "region_id")
    }


    dt <- dt %>%
        dplyr::select(-c(region_id, feature_id, overlap_w))

    cnames <- c(str_c(colnames(dt), ".A"), str_c(colnames(dt), ".B"))

    dt <- cbind(dt[gi@anchor1, ], dt[gi@anchor2, ])
    colnames(dt) <- cnames
    dt <- dt %>% dplyr::select(order(colnames(dt)))

    if (save_ambig) {
        dt <- dt %>%
            relocate(ambig.A, ambig.B, ambig_list.A, ambig_list.B, .after = last_col())
    } else {
        dt <- dt %>% dplyr::select(-c(ambig.A, ambig.B, ambig_list.A, ambig_list.B))
    }

    dt <- dt %>% data.frame()
    mcols(gi) <- cbind(mcols(gi)[orininal_columns], dt)
    return(gi)
}


#' Annotate RNA-RNA interactions as cis- and trans-
#'
#' @description
#' Annotated each entry gi object as cis, if the .A and .B arms correspond to
#' the same feature (i.e transcript_id or gene_id)
#' If the values are are equal, then annotated with value: `cis` = 1, If not equal
#' or `NA`: `cis` = 0
#'
#' @param gi `GInteractions` object containing two metadata columns as feature annotation
#' @param id_col_base base name of the feature id columns to use. Function will
#' look for <id_col_base>.A and <id_col_base>.B columns and compare them
#' @keywords internal
#' @return gi `GInteractions` object containing `cis` field with 0/1 values
.annotateCisTrans <- function(gi, id_col_base = "gene_id") {
    gi$cis <- 0

    name_col1 <- paste0(id_col_base, ".A")
    name_col2 <- paste0(id_col_base, ".B")

    gi$c1 <- mcols(gi)[, name_col1]
    gi$c2 <- mcols(gi)[, name_col2]


    if (length(gi[!is.na(gi$c1) & !is.na(gi$c2) &
        (gi$c1 == gi$c2)]) != 0) {
        gi[!is.na(gi$c1) & !is.na(gi$c2) &
            (gi$c1 == gi$c2)]$cis <- 1
        gi$c1 <- NULL
        gi$c2 <- NULL
        return(gi)
    } else {
        message("Could not do cis/trans annotation because columns content is missing ")
        return(gi)
    }
}

#' Helper function to add count data to metadata of `GInteractions`
#'
#' Merges the count dataframe and interactions metadata by `id_col`
#' If key is not found, in metadata throws error

#' @param gi `GInteractions`
#' @param df_counts dataframe with read counts
#' @param id_col key to use in merge
#' @return `GInteractions` with added counts
.addGeneCounts <- function(gi, df_counts, id_col = "gene_id") {
    df_counts <- df_counts %>% as.data.frame()
    df_all_cts <- as_tibble(data.frame("RNA" = unname(df_counts[1]), "n" = unname(df_counts[2])))

    id_columns <- str_c(id_col, c(".A", ".B"))

    if (!all(id_columns %in% colnames(mcols(gi)))) {
        stop(id_columns, " are not found in the input gi")
    }
    colnames(df_all_cts) <- c(id_columns[1], "gene_count.A")
    gi$gene_count.A <- left_join(
        as_tibble(mcols(gi)[id_columns[1]]),
        df_all_cts
    ) %>% pull(gene_count.A)
    colnames(df_all_cts) <- c(id_columns[2], "gene_count.B")
    gi$gene_count.B <- left_join(
        as_tibble(mcols(gi)[id_columns[2]]),
        df_all_cts
    ) %>% pull(gene_count.B)
    return(gi)
}


calculateLigationPvalues <- function(gi, df_counts, id_col = "gene_id") {
    df_counts <- df_counts %>% as.data.frame()
    df_all_cts <- as_tibble(data.frame("RNA" = unname(df_counts[1]), "n" = unname(df_counts[2])))

    id_columns <- str_c(id_col, c(".A", ".B"))

    if (!all(id_columns %in% colnames(mcols(gi)))) {
        stop(id_columns, " are not found in the input gi")
    }

    chim_df <- as_tibble(mcols(gi)) %>% dplyr::select(c(n_reads, !!!id_columns))
    colnames(chim_df) <- c("n", "Araw", "Braw")
    chim_df$chim_id <- seq_len(nrow(chim_df))

    unique_ids <- unique(df_all_cts$RNA)
    # filter not annotated chimeras and those for which we don't have record in counts df

    chim_df$AB <- apply(chim_df[c("Araw", "Braw")], 1, function(row) {
        paste(sort(row), collapse = "<>")
    })
    chimdf_save <- chim_df[c("chim_id", "AB")]

    chim_df <- chim_df %>%
        group_by(AB) %>%
        summarise(n_reads_chim_total = sum(n)) %>%
        ungroup()
    ssplit <- str_split_fixed(chim_df$AB, "<>", n = 2)
    chim_df$A <- ssplit[, 1]
    chim_df$B <- ssplit[, 2]
    chim_df <- chim_df %>%
        left_join(tibble("idA" = df_all_cts$RNA, "nA" = df_all_cts$n),
            by = c("A" = "idA")
        )
    chim_df <- chim_df %>%
        left_join(tibble("idB" = df_all_cts$RNA, "nB" = df_all_cts$n),
            by = c("B" = "idB")
        )
    chim_df <- chim_df %>% dplyr::filter((A %in% unique_ids) &
        (B %in% unique_ids))
    N_total <- sum(df_all_cts$n) + sum(chim_df$n_reads_chim_total)
    N_chim_total <- sum(chim_df$n_reads_chim_total)
    chim_df <- chim_df %>% mutate(
        Pa = nA / N_total,
        Pb = nB / N_total,
        Pab = ifelse(A != B, 2 * Pa * Pb, Pa * Pb)
    )
    chim_df$Pab_norm <- scale(chim_df$Pab, center = FALSE, scale = sum(chim_df$Pab))[, 1]
    chim_df$pval <- pbinom(chim_df$n_reads_chim_total, size = N_total, prob = chim_df$Pab_norm)
    chim_df$p.adj <- p.adjust(chim_df$pval, method = "BH")
    chimdf_save <- chimdf_save %>% left_join(chim_df[c("p.adj", "AB")], by = "AB")


    gi$p_val <- tibble("chim_id" = seq_len(length(gi))) %>%
        left_join(chimdf_save, by = "chim_id") %>%
        pull(p.adj)
    gi <- .addGeneCounts(gi, df_counts)

    return(gi)
}



.checkRNAduplexinstalled <- function() {
    returncode <- system2("RNAduplex", args = c("-h"), stdout = NULL, stderr = NULL)
    if (returncode != 0) {
        return(1)
    } else {
        return(0)
    }
}

.runRNAduplex <- function(RNA1, RNA2) {
    input <- paste(RNA1, RNA2, sep = "\n")
    result <- system2("RNAduplex", input = input, stdout = TRUE)
    return(result)
}

.getDuplexString <- function(gi, fafile) {
    s <- Biostrings::readBStringSet(fafile)
    s <- Biostrings::RNAStringSet(s)
    newnames <- str_squish(str_sub(names(s), 1, 5))
    names(s) <- newnames
    gi$seq1 <- s[get_arm_a(gi)]
    gi$seq2 <- s[get_arm_b(gi)]
    return(gi)
}


.getGCContent <- function(seqrna) {
    # seqrna is a RNAstringset based on the arm
    gc_freq <- Biostrings::letterFrequency(seqrna, c("G", "C"))
    gc_content <- rowSums(gc_freq) / width(seqrna) * 100
    return(gc_content)
}
