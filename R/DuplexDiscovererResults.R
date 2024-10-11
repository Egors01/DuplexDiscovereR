#' DuplexDiscovererResults
#'
#' A helper S4 class to store the results of the full
#' DuplexDiscovereR analysis.
#' This class contains the following output:
#' - `duplex_groups`:  clustered duplex groups.
#' - `chimeric_reads`: individual two-regions chimeric reads.
#'    Contains both clustered and unclustered reads.
#'    Clustered reads are linked to the duplex groups though 'dg_id' field in metadata
#' - `reads_classes`: dataframe parallel to the the input containing classification result and detected mapping type for each entry in the input
#' - `chimeric_reads_stats`: dataframe containing read type classification statistics
#' - `run_stats`: data frame containing statistics about the time and memory used by the pipeline
#' @details
#' Each output type has a corresponding accessor: 
#'  - [DuplexDiscovereR::dd_get_duplex_groups()]
#'  - [DuplexDiscovereR::dd_get_chimeric_reads()]
#'  - [DuplexDiscovereR::dd_get_reads_classes()]
#'  - [DuplexDiscovereR::dd_get_chimeric_reads_stats()]
#'  - [DuplexDiscovereR::dd_get_run_stats()]
#' @seealso 
#'  [DuplexDiscovereR::dd_get_duplex_groups()] ,
#'  [DuplexDiscovereR::dd_get_chimeric_reads()] ,
#'  [DuplexDiscovereR::dd_get_reads_classes()] ,
#'  [DuplexDiscovereR::dd_get_chimeric_reads_stats()] ,
#'  [DuplexDiscovereR::dd_get_run_stats()]
#'  
#' @slot duplex_groups \pkg{GInteractions} object with duplex groups
#' @slot chimeric_reads \pkg{GInteractions} object with chimeric reads
#' @slot reads_classes `tibble` (tbl_df) with read classification data.
#' @slot chimeric_reads_stats `tibble` (tbl_df)  read type statistics.
#' @slot run_stats `tibble` (tbl_df)  runtime and memory info
#' @rdname DuplexDiscovererResults
#' @return A `DuplexDiscovererResults` object. 
setClass("DuplexDiscovererResults",
    slots = c(
        duplex_groups = "GInteractions",
        chimeric_reads = "GInteractions",
        reads_classes = "tbl_df",
        chimeric_reads_stats = "tbl_df",
        run_stats = "tbl_df"
    )
)
#' Constructor for DuplexDiscovererResults class
#'
#' @param duplex_groups \pkg{GInteractions} object with duplex groups
#' @param chimeric_reads \pkg{GInteractions} object with chimeric reads
#' @param reads_classes `tibble` (tbl_df) with read classification data.
#' @param chimeric_reads_stats `tibble` (tbl_df)  read type statistics.
#' @param run_stats `tibble` (tbl_df)  runtime and memory info
#' @rdname DuplexDiscovererResults
#' @inherit dd_get_duplex_groups examples
#' @export
DuplexDiscovererResults <- function(duplex_groups,
    chimeric_reads,
    reads_classes,
    chimeric_reads_stats,
    run_stats) {
    new("DuplexDiscovererResults",
        duplex_groups = duplex_groups,
        chimeric_reads = chimeric_reads,
        reads_classes = reads_classes,
        chimeric_reads_stats = chimeric_reads_stats,
        run_stats = run_stats
    )
}

#' Show Method for DuplexDiscovererResults class
#'
#' This method provides a summary of the DuplexDiscovererResults object.
#' It prints `chimeric_reads_stats` followed by the `run_stats`.
#'
#' @param object A `DuplexDiscovererResults` object.
#' @return None. Prints a formatted summary.
#' @export
setMethod("show", "DuplexDiscovererResults", function(object) {
    df <- object@chimeric_reads_stats
    sname <- paste0("Sample name: ", df$sample_name[1])
    df$sample_name <- NULL
    transposed <- as.data.frame(t(df))
    transposed$feature <- rownames(transposed)
    colnames(transposed) <- c("count", "feature")
    transposed$count <- round(transposed$count, 2)
    transposed <- transposed[, c("feature", "count")]

    # Pretty print the transposed tibble
    message("DuplexDiscovereR Results")
    message(sname)
    message("Chimeric reads statistics:")
    print(transposed, row.names = FALSE)
    message("Use class accessors to retrieve output. I.e duplex_groups(resultsObject) ")
})

#' Accessor for `duplex_groups` slot
#'
#' Retrieves the value of the `duplex_groups` slot in a `DuplexDiscovererResults` object.
#' @param object A `DuplexDiscovererResults` object.
#' @return GInteractions object from the `duplex_groups` slot.
#' @docType methods
#' @rdname dd_get_duplex_groups
#' @export
#' @examples
#' 
#' # load example input 
#' data("RNADuplexesSmallGI")
#' data("RNADuplexesSampleData")
#' # run whole pipeline
#' result <- runDuplexDiscoverer(
#' data = SampleSmallGI,
#' junctions_gr = SampleSpliceJncGR,
#' anno_gr = SampleGeneAnnoGR,
#' sample_name = "run_example",
#' lib_type = "SE",
#' table_type = "STAR"
#' )
#' # access results 
#' show(result)
#' gi_clusters <- dd_get_duplex_groups(result)
#' gi_reads <- dd_get_chimeric_reads(result)
#' df_reads <- dd_get_reads_classes(result)
#' dd_get_reads_classes(result)
#' dd_get_run_stats(result)
setGeneric("dd_get_duplex_groups", function(object) standardGeneric("dd_get_duplex_groups"))

#' @rdname dd_get_duplex_groups
#' @docType methods
setMethod("dd_get_duplex_groups", "DuplexDiscovererResults", function(object) object@duplex_groups)


#' Accessor for `chimeric_reads` Slot
#'
#' Retrieves the value of the `chimeric_reads` slot in a `DuplexDiscovererResults` object.
#' @param object A `DuplexDiscovererResults` object.
#' @return GInteractions object from the `chimeric_reads` slot.
#' @docType methods
#' @rdname dd_get_chimeric_reads
#' @inherit dd_get_duplex_groups examples
#' @export
setGeneric("dd_get_chimeric_reads", function(object) standardGeneric("dd_get_chimeric_reads"))
#' @docType methods
#' @rdname dd_get_chimeric_reads
setMethod("dd_get_chimeric_reads", "DuplexDiscovererResults", function(object) object@chimeric_reads)

#' Accessor for `reads_classes` Slot
#'
#' Retrieves the value of the `reads_classes` slot in a `DuplexDiscovererResults` object.
#' @param object A `DuplexDiscovererResults` object.
#' @return `tibble` from the `reads_classes` slot.
#' @docType methods
#' @rdname dd_get_reads_classes
#' @inherit dd_get_duplex_groups examples
#' @export
setGeneric("dd_get_reads_classes", function(object) standardGeneric("dd_get_reads_classes"))

#' @docType methods
#' @rdname dd_get_reads_classes
setMethod("dd_get_reads_classes", "DuplexDiscovererResults", function(object) object@reads_classes)

#' Accessor for `chimeric_reads_stats` Slot
#'
#' Retrieves the value of the `chimeric_reads_stats` slot in a `DuplexDiscovererResults` object.
#' @param object A `DuplexDiscovererResults` object.
#' @return `tibble` from the `chimeric_reads_stats` slot.
#' @docType methods
#' @rdname dd_get_chimeric_reads_stats
#' @inherit dd_get_duplex_groups examples
#' @export
setGeneric("dd_get_chimeric_reads_stats", function(object) standardGeneric("dd_get_chimeric_reads_stats"))

#' @docType methods
#' @rdname dd_get_chimeric_reads_stats
setMethod("dd_get_chimeric_reads_stats", "DuplexDiscovererResults", function(object) object@chimeric_reads_stats)


#' Accessor for `run_stats` Slot
#'
#' Retrieves the value of the `run_stats` slot in a `DuplexDiscovererResults` object.
#' @param object A `DuplexDiscovererResults` object.
#' @return `tibble` from the `run_stats` slot.
#' @docType methods
#' @rdname dd_get_run_stats
#' @inherit dd_get_duplex_groups examples
#' @export
setGeneric("dd_get_run_stats", function(object) standardGeneric("dd_get_run_stats"))

#' @docType methods
#' @rdname dd_get_run_stats
setMethod("dd_get_run_stats", "DuplexDiscovererResults", function(object) object@run_stats)
