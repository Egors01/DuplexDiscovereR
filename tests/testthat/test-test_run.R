test_that("Executing all steps in single-call works", {
    suppressMessages({
        data(RNADuplexesSampleData)
        set.seed(123)
        res <- runDuplexDiscoverer(
            data = RNADuplexesRawChimSTAR,
            junctions_gr = SampleSpliceJncGR,
            anno_gr = SampleGeneAnnoGR,
            df_counts = RNADuplexesGeneCounts,
            sample_name = "test clustering",
            fafile = NULL,
            collapse_n_inter = 3,
            lib_type = "SE",
            table_type = "STAR"
        )
        gi_clusters <- res$gi_clusters
        gi_reads <- res$gi_reads
        df_reads <- res$df_reads

        n_reads_clustered <- sum(gi_clusters$n_reads)
        n_reads_clustered_stats <- df_reads %>%
            dplyr::filter(!is.na(dg_id)) %>%
            pull(n_reads) %>%
            sum()
        # dplyr::filter(map_type == "2arm") %>%
        n_reads_unclustered <- df_reads %>%
            dplyr::filter(is.na(dg_id)) %>%
            distinct(readname) %>%
            nrow()
        n_reads_passed_to_clustering <- df_reads %>%
            dplyr::filter(map_type == "2arm") %>%
            nrow()
        n_reads_passed_to_clustering_nodg <- df_reads %>%
            dplyr::filter(is.na(dg_id), map_type == "2arm") %>%
            nrow()
        n_clusters <- length(gi_clusters)
    })
    expect_equal(n_reads_clustered, n_reads_clustered_stats, label = "read numbers match between stats and returned GI")
    expect_equal(n_reads_clustered + n_reads_passed_to_clustering_nodg, n_reads_passed_to_clustering, label = "no reads lost before clustering")
    expect_equal(n_reads_clustered, 300, label = "correct number of reads clustered")
    expect_equal(n_clusters, 79, label = "correct number of clusters found")
})
