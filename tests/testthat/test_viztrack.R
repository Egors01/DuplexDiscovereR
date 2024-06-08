# generate example duplexes to plot
anchor1 <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(
        start = c(100, 600, 1100, 1600, 2100, 150, 400),
        end = c(200, 700, 1200, 1700, 2200, 250, 500)
    ),
    strand = "+"
)
anchor2 <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(
        start = c(300, 800, 1300, 1800, 2300, 1500, 1700),
        end = c(400, 900, 1400, 1900, 2400, 1600, 1800)
    ),
    strand = "+"
)

interactions <- GInteractions(anchor1, anchor2, mode = "strict")

# define plotting range, assigns some annotaiton
gr_region <- range(anchor1, anchor2)
interactions$anno_A <- sample(LETTERS, length(interactions))
interactions$anno_B <- interactions$anno_A

# create track

# add interactions which are not fully in plot range: outside the range
# or on different chromosome()
# one left (A) interaction arm outside of the plot,
# other on different chromosome
new_anchor1 <- GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges(
        start = c(10, 600),
        end = c(90, 700)
    ),
    strand = "+"
)
new_anchor2 <- GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges(
        start = c(1500, 1000),
        end = c(1600, 1200)
    ),
    strand = "+"
)

new_interactions <- GInteractions(new_anchor1, new_anchor2)
new_interactions$anno_A <- c("A.out", "A.out_chr")
new_interactions$anno_B <- c("B.in", "B.in")
all_interactions <- c(interactions, new_interactions)



test_that("Creating visualization track workds", {
    suppressMessages({
        track_a <- DuplexTrack(interactions, gr_region = gr_region, stacking = "dense")
        track_b <- DuplexDiscovereR::DuplexTrack(all_interactions,
            gr_region = gr_region
        )
    })
    expect_equal(length(range(track_a)), 14,
        label =
            "number or ranges on plot is correct:1"
    )
    expect_equal(length(range(track_b)), 16,
        label =
            "number or ranges on plot is correct:2"
    )
})
# for trackB check in-out-range calculation
expacted_in_names <- c(
    "1.A", "1.B", "2.A", "2.B", "5.A", "5.B", "8.B", "9.A", "6.A",
    "6.B", "7.A", "7.B", "3.A", "3.B", "4.A", "4.B"
)
expected_in_labels <- c(
    "1", "1", "1", "1", "1", "1", "1", "0", "1", "1", "1",
    "1", "1", "1", "1", "1"
)
names(expected_in_labels) <- expacted_in_names

# for trackA check plot in 'dense' to indirectly test AnotationTrack option worked
expected_rowids <- rep("1", 14)
names(expected_rowids) <- paste0(rep(seq(1:7), each = 2), ".", rep(c("A", "B"), 7))


test_that("Drawing visualization track workds", {
    set.seed(123)
    suppressMessages({
        track_a <- DuplexTrack(interactions,
            gr_region = gr_region,
            stacking = "dense"
        )
        track_b <- DuplexDiscovereR::DuplexTrack(all_interactions,
            gr_region = gr_region
        )
        r_a <- plotTracks(track_a, stacking = "dense")
        r_b <- plotTracks(track_b, stacking = "squish")
    })
    in_range <- tags(imageMap(r_b$DuplexTrack))$in_range
    rowids <- tags(imageMap(r_a$DuplexTrack))$row_plot_id


    expect_equal(in_range, expected_in_labels)
    expect_equal(rowids, expected_rowids)
})
