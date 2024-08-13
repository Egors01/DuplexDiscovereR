#' Run prediciton of RNA hybridization
#'
#' Calls RNAduplex from ViennaRNA to find base-pairs for every entry in the
#' input, throws a message and system warning if it is not installed
#'
#' @param gi `Ginteraction` with pairs of regions
#' @param fafile path to the .fasta file with genome
#'
#' @returns object parallel to input with added energy GC content,
#' dot-format base-pairings and lenghts of RNA hybrids
#' will return the input, if RNAhybrids cannot be run
#'
#' @export
#' @examples
#' sequence <- paste0(
#'     "AGCUAGCGAUAGCUAGCAUCGUAGCAUCGAUCGUAAGCUAGCUAGCUAGCAUCGAUCGUAGCUAGCAUCGAU",
#'     "CGUAGCAUCGUAGCUAGCUAGCUAUGCGAUU"
#' )
#'
#' # Save the sequence to a temp fasta file
#' fasta_file <- tempfile(fileext = ".fa")
#' chrom <- "test_chrA"
#' writeLines(c(">test_chrA", sequence), con = fasta_file)
#'
#' # Create the GInteraction object
#' # Define start and end positions for the base-pairing regions
#' regions <- data.frame(
#'     start1 = c(1, 11, 21, 31, 41),
#'     end1 = c(10, 20, 30, 40, 50),
#'     start2 = c(91, 81, 71, 61, 51),
#'     end2 = c(100, 90, 80, 70, 60)
#' )
#' # GRanges objects for the anchors
#' anchor1 <- GRanges(seqnames = chrom, ranges = IRanges(start = regions$start1, end = regions$end1))
#' anchor2 <- GRanges(seqnames = chrom, ranges = IRanges(start = regions$start2, end = regions$end2))
#' interaction <- GInteractions(anchor1, anchor2)
#' # predict hybrids
#' # In case ViennaRNA is installed
#' \dontrun{
#' getRNAHybrids(interaction, fasta_file)
#' }
getRNAHybrids <- function(gi, fafile) {
    # exit if not installed
    if (.checkRNAduplexinstalled() != 0) {
        return(gi)
    }

    sq <- Biostrings::readBStringSet(fafile)
    sq <- Biostrings::DNAStringSet(sq)
    sq <- Biostrings::RNAStringSet(sq)
    if (!all(as.character(seqnames(regions(refresh_gi(gi)))) %in% names(sq))) {
        stop("Calling hybrids: Seqnames in fasta file do not contain all seqnames of the input gi ")
    }

    gi$GC_content.A <- .getGCContent(sq[get_arm_a(gi)])
    gi$GC_content.B <- .getGCContent(sq[get_arm_b(gi)])

    gr1 <- get_arm_a(gi)
    gr2 <- get_arm_b(gi)
    seq1 <- sq[gr1]
    seq2 <- sq[gr2]

    rr <- mapply(.runRNAduplex, seq1, seq2)

    dots <- str_split_i(rr, " ", 1)
    energies <- as.numeric(sub(".*\\((-?[0-9]+\\.[0-9]+)\\).*", "\\1", rr))

    matches <- str_match(rr, "\\s(\\d+),(\\d+)\\s+:\\s+(\\d+),(\\d+)\\s")
    len1 <- as.numeric(matches[, 3]) - as.numeric(matches[, 2])
    len2 <- as.numeric(matches[, 5]) - as.numeric(matches[, 4])
    gi$E <- energies
    gi$hyb_len_1 <- len1
    gi$hyb_len_2 <- len2
    gi$dots <- dots
    return(gi)
}
