n_dgs_expected <- 1830
chrlen_ln <- 50818468
test_that("write to sam file works", {
    suppressMessages({
        data("RNADuplexesSampleData")
        length(RNADuplexSampleGI)
        table(is.na(RNADuplexSampleGI$dg_id))
        tmpf <- tempfile(".sam")
        writeGiToSAMfile(
            gi_coords = RNADuplexSampleGI,
            id_column = "dg_id",
            file_out = tmpf,
            distance_chim_junction = 1e5,
            genome = "hg38"
        )
        sam_body <- read.delim(tmpf, header = FALSE, comment.char = "@")
        dgs_assigned <- as.integer(stringr::str_split_i(sam_body$V13, "DG:Z:", 2))
        dgs_count_tab <- table(dgs_assigned)


        sam_header <- read.delim(tmpf, header = FALSE, nrows = 4)
        chrlen_real <- gsub(".*LN:(\\d+).*", "\\1", paste0(sam_header$V3, collapse = ""))
        chrlen_real <- as.integer(chrlen_real)
    })
    expect_equal(dgs_count_tab[["0"]], n_dgs_expected, label = "NA dgs assigned to zero")
    expect_equal(chrlen_real, chrlen_ln, label = "Chrlen is written")
})
