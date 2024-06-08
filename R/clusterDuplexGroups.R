#' Cluster RNA duplexes in `GInteractions` object
#'
#' @description
#' Main method to find duplex groups from the individual interactions
#'
#' @details
#' Accepts or creates the connections graphdf dataframe, creates graph with
#' `igraph` package, uses community detection algoritm to call clusters.
#' New field `dg_id` is added to label the clusters (duplex groups).
#' If no community is found for the read, `dg_id` is `NA`
#'
#' @param gi `GInteractions` object
#' @param min_arm_ratio For graph creation only. Span-to-overlap ratio threshold. If smaller than this value, then edge is not drawn
#' @param maxgap For graph creation only. Max shift between arms starts and ends for pair of overlapping reads
#' @param minovl For graph creation only. Minimum required overlap between either arm for pair of overlapping reads
#' Other optional arguments, which are not relevant, unless user want to modify clustering
#' weights or modify clustering in some other way
#' @param graphdf Optional. Dataframe representing connection edges between entries in gi
#' If not provided, graphdf is created inside the function
#' @param id_column Optional. Column name in the  `GInteractions` metadata,
#'  which was used to index temporary duplex groups, if they are present
#' @param weight_column  Optional. If graphdf is provided,  field to use for weight overlaps
#' @param fast_greedy  Optional. Run the fast_greedy algorithm instead of Louvain.
#' Can speed up calcualtion for the large graphs.
#' @param decompose Decompose graph into separate sub-graphs before clustering.
#  Increases runtime, sometimes results in higher number of overlapping DG.
#' @param id_columns_grapdf Column in the graph dataframe, which was used for index
#' @param dump_graph For debug. Export the graph elements. not used
#' @param dump_path For debug. PArt to export the graph elements. not used
#'
#' @return `GInteractions` object with new \code{dg_id} column
#' @export
#' @examples
#' data("RNADuplexesSampleData")
#' # run preprocessing and filtering
#' preproc_df <- runDuplexDiscoPreproc(RNADuplexesRawBed, table_type = "bedpe")
#' preproc_gi <- makeGiFromDf(preproc_df)
#' preproc_gi <- classifyTwoArmChimeras(preproc_gi,
#'     min_junction_len = 5,
#'     junctions_gr = SampleSpliceJncGR, max_sj_shift = 10
#' )
#' # collapse duplicates
#' gi <- collapseIdenticalReads(preproc_gi)$gi
#' # run global clustering
#' gi <- clusterDuplexGroups(gi)
#' # check dg_ids
#' table(is.na(gi$dg_id))
clusterDuplexGroups <- function(gi, graphdf = NULL, maxgap = 40,
    minovl = 10,
    id_column = "duplex_id",
    weight_column = "weight",
    fast_greedy = FALSE,
    decompose = FALSE,
    id_columns_grapdf = paste(id_column, c(1, 2), sep = "."),
    min_arm_ratio = 0.3,
    dump_graph = FALSE,
    dump_path = "") {
    gi$dg_id <- NULL

    if (is.null(graphdf)) {
        message("Computing overlaps on-the-fly")
        gi$idx_tmp <- seq_len(length(gi))
        id_column <- "idx_tmp"
        graphdf <- computeGISelfOverlaps(gi,
            id_column = id_column, maxgap = maxgap,
            minovl = minovl
        )
        if (nrow(graphdf) == 0) {
            message("No overlap found. No DG found in input object")
            gi$dg_id <- NA
            return(gi)
        }
        id_columns_grapdf <- paste(id_column, c(1, 2), sep = ".")
        weight_column <- "weight"
        graphdf <- graphdf %>% dplyr::filter(
            ratio.A >= min_arm_ratio,
            ratio.B >= min_arm_ratio
        )
    } else {
        # message("Climeric duplex id column: ",id_column)
        # message("Using connections dataframe with indexes: ",id_columns_grapdf[1],"-",id_columns_grapdf[2])
        # message("Weights in: ",weight_column)
        if (!(id_column %in% colnames(mcols(gi)))) {
            stop(id_column, "not found in provided GenomicInteractions onject")
        } else {
            if (any(is.na(mcols(gi)[id_column]))) {
                stop(paste0("Undefined values in index field found before duplex clustering. Index column used: ", id_column, collapse = ""))
            }
        }
    }
    graphdf$weight <- graphdf[weight_column] %>% pull()
    graphdf$vert_id1 <- graphdf[id_columns_grapdf[1]] %>% pull()
    graphdf$vert_id2 <- graphdf[id_columns_grapdf[2]] %>% pull()
    graphdf <- graphdf %>% dplyr::select(vert_id1, vert_id2, weight)
    graph_vertices <- unique(c(graphdf$vert_id1, graphdf$vert_id2))

    g <- igraph::graph_from_data_frame(graphdf, directed = FALSE, vertices = graph_vertices)
    igraph::edge.attributes(g)$weight <- round(graphdf$weight, 3)
    g <- igraph::simplify(g, edge.attr.comb = list(weight = "sum"), remove.multiple = TRUE)


    if (decompose) {
        # here we split into subgraphs and cluster each
        comps <- igraph::decompose(g)
        message("Whole-transciptome graph was split into ", length(comps), " components")
        pb <- txtProgressBar(min = 1, max = length(comps), style = 3)
        for (i in seq_along(comps)) { # For each subgraph comps[[i]]
            wt <- igraph::cluster_louvain(comps[[i]], weights = E(comps[[i]])$weight)
            igraph::V(comps[[i]])$cluster_group <- str_c(i, "_", igraph::membership(wt))
            setTxtProgressBar(pb, i)
        }
        close(pb)
        df_clustered <-
            tibble(
                dg_id_raw = unlist(lapply(comps, function(x) {
                    V(x)$cluster_group
                })),
                vert_id = as.integer(unlist(lapply(comps, function(x) {
                    names(V(x))
                })))
            ) %>%
            mutate(idx = as.integer(vert_id)) %>%
            group_by(dg_id_raw) %>%
            mutate(
                dg_id = cur_group_id(),
                n_el = n()
            ) %>%
            ungroup() %>%
            mutate(cmp = str_split_i(dg_id_raw, "_", 1)) %>%
            group_by(cmp) %>%
            mutate(n_clust = n_distinct(dg_id)) %>%
            ungroup()
        if (dump_graph) {
            grres <- list()
            grres$comps <- comps
            grres$df_clustered <- df_clustered
            save(grres, file = dump_path)
        }
    } else {
        # here we do clustering genome-wide
        if (!fast_greedy) {
            cluster_lou <- igraph::cluster_louvain(g)
        } else {
            cluster_lou <- igraph::cluster_fast_greedy(g)
        }
        igraph::V(g)$cluster_group <- igraph::membership(cluster_lou)

        # create output df for clustered duplexes
        df_clustered <- tibble(
            dg_id = as.integer(igraph::V(g)$cluster_group),
            vert_id = as.integer(names(igraph::V(g)))
        ) %>%
            arrange(dg_id)
        remove(g)
        remove(graphdf)
        gc()
    }

    if ("idx_tmp" %in% colnames(mcols(gi))) {
        gi$idx <- NULL
    }

    gi$dg_id <- left_join(
        tibble(vert_id = as.integer(mcols(gi)[, id_column])),
        df_clustered[, c("vert_id", "dg_id")],
        by = "vert_id"
    ) %>%
        pull(dg_id)

    return(gi)
}
