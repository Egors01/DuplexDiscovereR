test_that("Load test data object run small clustering works", {
    suppressMessages({
        data("RNADuplexesSmallGI")
        SampleSmallGI$n_reads <- 1
        gi <- clusterDuplexGroups(SampleSmallGI, maxgap = 40, minoverlap = 10, min_arm_ratio = 0.3)
    })
    expect_equal(length(collapse_duplex_groups(gi)), 3)
    expect_equal(length(collapse_duplex_groups(gi, return_unclustered = T)), 6)
})
test_that("Load STAR format and classification clustering works", {
    suppressMessages({
        set.seed(123)
        test_reads_star <- read.table(system.file("extdata", "test_SPLASH_Chimeric.out.junction", package = "DuplexDiscoverer"), header = T)
        df_chim <- runDuplexDiscoPreproc(data = test_reads_star, library_type = "SE", table_type = "STAR", keep_metadata = F)
        typetable <- table(df_chim$map_type)
        two_arm_N <- typetable[["2arm"]]
        mm_N <- typetable[["multi_map"]]
        ms_N <- typetable[["multi_split"]]
    })
    expect_equal(two_arm_N, "expected" = 2090, label = "mapping type classification works 1 ")
    expect_equal(mm_N, "expected" = 218, label = "mapping type classification works 2")
    expect_equal(ms_N, "expected" = 285, label = "mapping type classification works 3")
})
test_that("Load bedpe format and crude clustering works", {
    suppressMessages({
        df_bedpe <- read.table(system.file("extdata", "test_SPLASH_DuplexesRaw.bedpe", package = "DuplexDiscoverer"))
        colnames(df_bedpe) <- get_colnames_and_types_for_input("BEDPE_COLNAMES")
        set.seed(123)
        df_chim <- runDuplexDiscoPreproc(data = df_bedpe, library_type = "SE", table_type = "bedpe", keep_metadata = F)
        df_chim <- df_chim %>% dplyr::filter(map_type == "2arm")
        gi_clusters <- collapse_duplex_groups(clusterDuplexGroups(makeGiFromDf(df_chim)))
        n_reads_clustered <- sum(gi_clusters$n_reads)
    })
    expect_equal(n_reads_clustered, 1162, label = "clustering bedpe produced expected results")
})
test_that("Load GInteractions format and crude clustering works", {
    suppressMessages({
        data("RNADuplexesSampleData")
        set.seed(123)
        df_chim <- runDuplexDiscoPreproc(data = RNADuplexSampleGI, keep_metadata = F)
        df_chim <- df_chim %>% dplyr::filter(map_type == "2arm")
        gi_clusters <- collapse_duplex_groups(clusterDuplexGroups(makeGiFromDf(df_chim)))
        n_reads_clustered <- sum(gi_clusters$n_reads)
    })
    expect_equal(n_reads_clustered, 1162, label = "clustering bedpe produced expected results")
})
