#' Executes all steps of DuplexDiscovereR pipeline
#' @description
#' Generates GInteractions object with duplex groups from the STAR
#' Chimeric.out.junction or bedpe file.
#' Classifies reads, annotates reads by overlap with the gene or transcript
#' features, calculates p-values and hybridization energies.
#' Additionally, returns mappings from duplex groupd back to genes.
#' @details
#' This is a main function to do the initial discovery of the RNA duplexes after
#' the chimeric read mapping. It wraps following procedures:
#' - Classifies the input reads by the mapping type.
#' Keeps 2-arm chimeric reads for downstream analysis
#' - Compares 2arm duplex reads against provided splice junctions
#' - Classifies 2arm duplexes into spurious self-overlapping, splice junction categoris
#' - Performs clustering of the remaining reads into duplex groups
#'    - Collapses identically mapped reads
#'    - Collapses closely located reads, almost identical reads
#'    - Finds duplex groups throughout whole  data set
#' - Annotates duplex groups with genomic features if annotation is provided
#' - Calculates p-values if gene counts and annotation are provided
#' - Calculates hybridization energies if path to the .fasta file is provided
#' @param data dataframe-like object with the split reads. Output of Chimeric.out.junction or dataframe with fileds defined by bedpe format:
#' c("chromA","startA",'endA',"chromB",'startB','endB','readname','flag','strandA','strandB', ... )
#' Alternatively, `GInteractions` object
#' @param table_type one in c("STAR","bedpe") Defines the type of the input dataframe.
#' ignored if input data is `GInteractions`
#' @param junctions_gr \pkg{GRanges} object with the splice junction coordinates
#' @param anno_gr \pkg{GRanges} object to use for the annotation of the interactions.
#' The c('gene_id','gene_name','gene_types') columns in anno_gr are used by default. Optional
#' @param fafile path to the genome .fasta file. Used to calculate hybridization energy with *RNADuplex*. Sequence names should correspond to the sequences from which the mapping index was created. Optional
#' @param df_counts A two- column dataframe with counts to use for p-value calculation. The first column should match the 'gene_id' feature in anno_gr. The second column is the respective count. Optional
#' @param sample_name A name of the sample, used for assembling the analysis statistics dataframe
#' @param lib_type one in c('SE','PE'). Type of the seqeuncing library. Default is 'SE'
#' @param max_sj_shift Maximum shift between either donor and acceptor splice sites and chimeric junction coordinates to count chimeric junction as splice junction
#' @param min_junction_len a minimum allowed distance between chimeric arms for the read input.
#' Reads with the junction closer than \code{min_junction_len} are annotated as '2arm_shot' and not clustered to duplex groups
#' @param max_gap Parameter for read clustering. Minimum required shift between start and end coordinates of arms for pair of overlapping chimeric reads.
#' If the shift is longer than \code{max_gap} for either arm, then total read overlap between those reads is zero.
#' @param min_overlap Parameter for read clustering. Minimum required overlap to for either arm (A or B) for pair of chimeric reads.
#' @param min_arm_ratio Parameter for read clustering.
#' If the overlap-to-span ratio for either arm (A or B) for pair of chimeric reads is less than \code{min_arm_ratio}, then the total overlap for this pair is set to zero.
#' @param collapse_n_inter  Parameter for read clustering  (iterative step). Number of iterations to repeat step of collapsing of the highly similar chimeric reads.
#' Increasing this from i.e 0 to 5 reduces clustering time and memory for the libraries with many overlapping reads.
#' @param gap_collapse_similar  Parameter for read clustering (iterative step). Analogous to the max_gap, but applied \code{collapse_n_inter} times during the iterative merging step.
#' Reduce this to 1 or 2 to lower RAM usage for clustering the library with many similar reads.
#' @return a list with the   following keys
#' \describe{
#'   \item{gi_clusters}{ \pkg{GInteractions} object with called duplex groups }
#'   \item{gi_reads}{ \pkg{GInteractions} object with non-collapsed reads, where 'dg_id' column attributes reads to the duplex groups}
#'   \item{df_reads}{ \pkg{tibble} parallel to the the input dataframe, annotated with read categories and duplex groups    }
#'   \item{df_stats}{ \pkg{tibble} dataframe with the infromation about the run}
#' }
#' @export
#'
#' @examples
#'
#' library(DuplexDiscovereR)
#' # load data
#' data("RNADuplexesSampleData")
#' res <- runDuplexDiscoverer(
#'     data = RNADuplexesRawChimSTAR,
#'     junctions_gr = SampleSpliceJncGR,
#'     anno_gr = SampleGeneAnnoGR,
#'     df_counts = RNADuplexesGeneCounts,
#'     sample_name = "test clustering",
#'     fafile = NULL,
#'     collapse_n_inter = 3,
#'     lib_type = "SE",
#'     table_type = "STAR"
#' )
#' names(res)
runDuplexDiscoverer <- function(data,
    table_type = "",
    junctions_gr = NULL,
    anno_gr = NULL,
    fafile = NULL,
    df_counts = NULL,
    sample_name = "sample",
    lib_type = "SE",
    min_junction_len = 5,
    max_gap = 50,
    min_arm_ratio = 0.1,
    min_overlap = 10,
    max_sj_shift = 10,
    gap_collapse_similar = 3,
    collapse_n_inter = 5) {
    memstart <- sum(data.frame(gc(reset = TRUE))[, 6])
    start_time <- Sys.time()
    # STEP 1 pre-process------
    time1 <- Sys.time()
    df <- runDuplexDiscoPreproc(data,
        table_type = table_type,
        library_type = lib_type,
        keep_metadata = TRUE
    )

    n_reads_initial <- sum(df$n_reads)
    message("Number of input alignments: ", n_reads_initial)
    time2 <- Sys.time()
    time_preproc <- round(as.numeric(difftime(time2, time1,
        units = "secs"
    )), 3)


    # STEP 2 select only 2-arm type chimeras, find SJ, find "too short" chimeras
    time1 <- Sys.time()
    # 2.1 filter out multi spli and multimap aln
    single_gap_df <- df %>%
        dplyr::filter(map_type == "2arm")
    # 2.1a create stats
    read_stats_df <- df %>%
        dplyr::select(readname, n_reads, map_type) %>%
        mutate(read_id = c(seq_len(nrow(df)))) %>%
        relocate(read_id, .after = n_reads)

    # 2.2 convert to GInteractions and mark short/overlapping/splice junction reads
    big_gi <- makeGiFromDf(single_gap_df)
    if (!is.null(junctions_gr)) {
        big_gi <- classifyTwoArmChimeras(
            gi = big_gi,
            min_junction_len = min_junction_len,
            junctions_gr = junctions_gr,
            max_sj_shift = max_sj_shift
        )
    } else {
        message(" Splice junction refrence is not provided")
        big_gi <- getChimericJunctionTypes(big_gi, normal_gap_threshold = min_junction_len)
        big_gi$splicejnc <- 0
    }

    # 2.2a save stats
    read_stats_df <- left_join(read_stats_df,
        as_tibble(mcols(big_gi)[c(
            "read_id", "junction_type",
            "splicejnc"
        )]),
        by = "read_id"
    ) %>%
        mutate(
            read_type = ifelse(is.na(junction_type),
                map_type, as.character(junction_type)
            ),
            read_type = ifelse(splicejnc == 0 | is.na(junction_type),
                read_type, "2arm_sj"
            )
        )
    time2 <- Sys.time()
    time_classify <- round(as.numeric(difftime(time2, time1,
        units = "secs"
    )), 3)

    # 2.3 keep only chimeric reads with 2 arms and not too short and not splice jnc
    big_gi$keep <- (big_gi$junction_type == "2arm") & (big_gi$splicejnc == 0)
    gi_2arm <- big_gi[big_gi$keep == TRUE]

    # STEP 3 clustering---------
    message("--- clustering ---")
    time1 <- Sys.time()

    # 3.1 prepare clustering: reduce complexity -----------
    message("--- collapsing identical reads ---")
    res_collapse_ident <- collapseIdenticalReads(gi_2arm)

    # 3.1a get results of collapse: get new gi object, update read stats
    gi <- res_collapse_ident$gi_collapsed
    read_stats_df <- left_join(read_stats_df,
        res_collapse_ident$stats_df,
        by = "read_id"
    ) %>%
        mutate(collapse_identical = 1)

    # 3.2 collapse (cluster) similar (less than 1..5 nt shift) reads
    # because thye will nodes be higly connected nodes"
    if (collapse_n_inter != 0) {
        message("--- iteratively collapse similar reads  --- ")
        message("minimum shift is  :", gap_collapse_similar, " nt")
        res <- collapseSimilarChimeras(gi, read_stats_df,
            maxgap = gap_collapse_similar,
            niter = collapse_n_inter,
            minoverlap = 20
        )
        # get results
        gi <- res$gi_updated
        read_stats_df <- res$stats_df
    }

    # 3.3 Clustering on the whole-genome -----
    message("--- calculating total read overlaps ---")

    graphdf_fast <- computeGISelfOverlaps(gi,
        maxgap = max_gap,
        id_column = "duplex_id",
        minoverlap = min_overlap
    )

    if ((nrow(graphdf_fast) == 0)) {
        message("Consider re-running pipeline with disabled extra collapse by
            setting collapse_n_inter = 0 ")
        stop("Global clustering cannot be called not run because there no
    further DG merge is possible")
    }

    # scale weights from (min:max): (0,2) to (0,1)
    graphdf_fast$weight <- scales::rescale(graphdf_fast$weight, to = c(0, 1))
    # prune weak edges
    graphdf_fast <- graphdf_fast %>% dplyr::filter(
        ratio.A >= min_arm_ratio,
        ratio.B >= min_arm_ratio
    )
    message("--- finding duplex groups  ---")
    # clustering will add one column "dg_id"
    gi_fast <- clusterDuplexGroups(gi, graphdf = graphdf_fast, decompose = FALSE)
    # add dg_id for duplexes which are aggregated locally, but not clustered globally
    gi_fast <- .addDGidsForTmpDGs(gi_fast)

    # use it to collapse duplexes into suplex groups
    gi_final <- collapse_duplex_groups(gi_fast,
        return_unclustered = FALSE,
        return_collapsed = TRUE,
        keep_meta = FALSE
    )
    # update read stats
    read_stats_df <- left_join(read_stats_df, .DGIdToDuplexId(gi_fast),
        by = "duplex_id"
    )

    dt_2arm <- left_join(tibble("read_id" = gi_2arm$read_id), read_stats_df,
        by = "read_id"
    ) %>%
        group_by(dg_id) %>%
        mutate("n_reads_dg" = n()) %>%
        ungroup() %>%
        mutate(n_reads_dg = ifelse(!is.na(dg_id), n_reads_dg, 0)) %>%
        dplyr::select(dg_id, duplex_id, n_reads_dg)

    gi_2arm$dg_id <- dt_2arm$dg_id
    gi_2arm$duplex_id <- dt_2arm$duplex_id
    gi_2arm$n_reads_dg <- dt_2arm$n_reads_dg
    gi_2arm$was_clustered <- 1

    # add the reads, which were dropped before clustering
    gi_2arm_full <- c(gi_2arm, big_gi[big_gi$keep == FALSE])
    gi_2arm_full$was_clustered <- ifelse(!is.na(gi_2arm_full$was_clustered), 1, 0)

    time2 <- Sys.time()
    time_clust <- round(as.numeric(difftime(time2, time1,
        units = "secs"
    )), 3)

    # STEP 4 Annotation...-------
    time1 <- Sys.time()
    if (!is.null(anno_gr)) {
        message("--- annotation --- ")
        gi_final <- annotateGI(gi_final, anno_gr)
        gi_final <- .annotateCisTrans(gi_final)
        not_annotated <- sum(as.integer(is.na(gi_final$gene_id.A) | is.na(gi_final$gene_id.B)))
        not_annotated_full <- sum(as.integer(is.na(gi_final$gene_id.A) & is.na(gi_final$gene_id.B)))
        annotated <- length(gi_final) - not_annotated
        message("N annotated duplex groups: ", annotated)
        message("N duplex groups with at least one arm missing annotaton: ", not_annotated)
        message("N duplex groups with at both arms missing annotaton: ", not_annotated_full)

        if (!is.null(df_counts)) {
            message("--- computing random ligation p-values ---")
            # gi_final = calculateLigationPvalues(gi_final,df_counts)
            gi_final <- calculateLigationPvalues(gi_final, df_counts)
        }
    } else {
        message("No annotation provided")
    }
    time2 <- Sys.time()
    time_anno <- round(as.numeric(difftime(time2, time1,
        units = "secs"
    )), 3)

    sumreads <- as_tibble(data.frame(mcols(gi_fast))) %>%
        dplyr::filter(!is.na(dg_id)) %>%
        group_by(dg_id) %>%
        summarise(nreads = sum(n_reads))

    widthA <- width(get_arm_a(gi_final))
    widthB <- width(get_arm_b(gi_final))
    meanA <- mean(widthA)
    medA <- median(widthA)
    meanB <- mean(widthB)
    medB <- median(widthB)

    end_time <- Sys.time()
    time_diff <- round(as.numeric(difftime(end_time, start_time,
        units = "secs"
    )), 3)
    memend <- sum(data.frame(gc())[, 6])
    memdiff <- round((memend - memstart) / 1024, 3)
    maxmem <- round((memend) / 1024, 3)
    mem_initial <- round((memstart) / 1024, 3)

    time1 <- Sys.time()
    if (!is.null(fafile)) {
        message("--- calculating hybridization ---")
        gi_final <- getRNAHybrids(gi_final, fafile)
    }
    time2 <- Sys.time()
    time_hyb <- round(as.numeric(difftime(time2, time1,
        units = "secs"
    )), 3)

    summary_stats_case <- tibble(
        "sample_name" = sample_name,
        "n_reads" = n_reads_initial,
        "exec_time" = time_diff,
        "memstart" = memstart,
        "memend" = memend,
        "memdiff" = round((memend - memstart) / 1024, 3),
        "maxmem" = round((memend) / 1024, 3),
        "mem_initial" = round((memstart) / 1024, 3),
        "time_preproc" = time_preproc,
        "time_classify" = time_classify,
        "time_clust" = time_clust,
        "time_hybrids" = time_hyb,
        "time_anno" = time_anno,
        "meanA" = meanA,
        "medianA" = medA,
        "meanB" = meanB,
        "medianB" = medB
    )

    gc(reset = TRUE)
    results <- list()
    results$gi_clusters <- gi_final
    results$gi_reads <- gi_2arm_full
    results$df_runstats <- summary_stats_case
    results$df_reads <- read_stats_df
    message("finished")
    return(results)
}
