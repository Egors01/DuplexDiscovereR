# Use artificial interactions
chrom <- "chr1"
start1 <- c(1, 5, 6, 5, 45, 45, 30, 30, 90, 90, 110, 10, 45, 48, 47)
end1 <- start1 + 10
start2 <- c(100, 100, 95, 100, 130, 130, 70, 70, 120, 120, 200, 200, 70, 68, 67)
end2 <- start2 + 10
expected_ovl_1stgroup <- 3
anchor1 <- GRanges(seqnames = chrom, ranges = IRanges(start = start1, end = end1))
anchor2 <- GRanges(seqnames = chrom, ranges = IRanges(start = start2, end = end2))

interaction <- GInteractions(anchor1, anchor2, mode = "strict")

n <- length(interaction)
group_size <- ceiling(n / 3) # should be 5
set.seed(123)
interaction$groups <- sample(rep(1:3, length.out = n))
# Split the dataframe into three groups
group1 <- interaction[interaction$groups == 1]
group2 <- interaction[interaction$groups == 2]
group3 <- interaction[interaction$groups == 3]
test_that("comparing samples works", {
    suppressMessages({
        a <- list("sample1" = group1, "sample2" = group2, "sample3" = group3)
        res <- compareMultipleInteractions(a,
                                           min_ratio = 0.0,
                                           minoverlap = 1)
        n_reads_real <- sum(res$gi_all$n_reads)
        len_superset <- length(res$gi_all)
        df_res <- res$dt_upset

        real_olv_1st <- sum(df_res[df_res$id == 1, c("sample1", "sample2", "sample3")])
        real_olv_six <- sum(df_res[df_res$id == 6, c("sample1", "sample2", "sample3")])
    })

    expect_equal(n_reads_real, 15, label = "superset size counted reads")
    expect_equal(len_superset, 7, label = "superset size is correct")
    expect_equal(real_olv_1st, expected_ovl_1stgroup,
        label = "test overlap correct:1"
    )
    expect_equal(real_olv_six, 1,
        label = "test overlap correct:6"
    )
})
