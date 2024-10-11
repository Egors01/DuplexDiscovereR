# Function to generate strings based on the 'length' column value
set_hg38_seqlengths <- function(gr) {
    # Define the sequence lengths for the hg38 assembly
    hg38_seqlengths <- c(
        chr1 = 248956422, chr2 = 242193529, chr3 = 198295559,
        chr4 = 190214555, chr5 = 181538259, chr6 = 170805979,
        chr7 = 159345973, chr8 = 145138636, chr9 = 138394717,
        chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
        chr13 = 114364328, chr14 = 107043718, chr15 = 101991189,
        chr16 = 90338345, chr17 = 83257441, chr18 = 80373285,
        chr19 = 58617616, chr20 = 64444167, chr21 = 46709983,
        chr22 = 50818468, chrX = 156040895, chrY = 57227415,
        chrM = 16569
    )

    # Intersect the provided GRanges object seqlevels with the hg38 seqlevels
    common_seqlevels <- intersect(seqlevels(gr), names(hg38_seqlengths))

    # Set the seqlengths for the GRanges object
    seqlengths(gr) <- hg38_seqlengths[common_seqlevels]

    return(gr)
}

set_hg19_seqlengths <- function(gr) {
    # Define the sequence lengths for the hg19 assembly
    hg19_seqlengths <- c(
        chr1 = 249250621, chr2 = 243199373, chr3 = 198022430,
        chr4 = 191154276, chr5 = 180915260, chr6 = 171115067,
        chr7 = 159138663, chr8 = 146364022, chr9 = 141213431,
        chr10 = 135534747, chr11 = 135006516, chr12 = 133851895,
        chr13 = 115169878, chr14 = 107349540, chr15 = 102531392,
        chr16 = 90354753, chr17 = 81195210, chr18 = 78077248,
        chr19 = 59128983, chr20 = 63025520, chr21 = 48129895,
        chr22 = 51304566, chrX = 155270560, chrY = 59373566,
        chrM = 16571
    )
    common_seqlevels <- intersect(seqlevels(gr), names(hg38_seqlengths))

    # Set the seqlengths for the GRanges object
    seqlengths(gr) <- hg38_seqlengths[common_seqlevels]

    return(gr)
}

#' Write reads to sam file
#' @description
#' Writes interactions to the sam file for visualization in extrnal browsers.
#' Takes input as GInteractions object containing reads or duplex groups.
#' @param gi_coords input Ginteraction object
#' @param genome DNAStringSet object with the reference genome
#' @param distance_chim_junction maximum distance between input
#' duplex groups/reads, which will be represented as the single-line in .sam file.
#' Junction will be output as N- gap. For the interactions with longer distances,
#' chimeric junction will be represented as MR:Z:i tag
#' @param read_name_column character field, pointing out to read names.
#' Read names are generated automatically if not provided.
#' @param id_column character name of the field containing integer duplex group ids.
#'  NA are replaced with zeros
#' @param genome character. Genome version. Required for the retrieval of sequence lengths
#' for sam file header-  SQ and SN tags.
#' For convenience, hg38 and hg19 chromosome lengths will be assigned
#'  automatically.
#' If the value is not in c('hg38','hg19'), seqlengths will be looked for be in
#' attribute in seqlengths() of regions(gi_coords)
#' @param file_out path to write output file
#' @returns no object is returned
#' @export
#' @examples
#' # Load test data
#' data("RNADuplexesSampleData")
#' # if the input is read-based, it should have integer duplex group ids
#' # here, we have 2090 reads
#' length(RNADuplexSampleGI)
#' # among them 300 reads does not belong to any DG
#' # missing ids will be converted to 0
#' table(is.na(RNADuplexSampleGI$dg_id))
#' tmpf <- tempfile(".sam")
#' writeGiToSAMfile(
#'     gi_coords = RNADuplexSampleGI,
#'     id_column = "dg_id",
#'     file_out = tmpf,
#'     distance_chim_junction = 1e5,
#'     genome = "hg38"
#' )
writeGiToSAMfile <- function(gi_coords, file_out,
    distance_chim_junction = 10000,
    read_name_column = "readname",
    id_column = "dg_id",
    genome = "") {
    # read names from metadata or generate
    if (read_name_column %in% colnames(mcols(gi_coords))) {
        gi_coords$rname <- mcols(gi_coords)[, read_name_column]
    } else {
        gi_coords$rname <- paste0("read_", seq_len(length(gi_coords)))
    }
    # id column - fill with 0 if there is a NA DG ids i
    if (id_column %in% colnames(mcols(gi_coords))) {
        gi_coords$GRP_id <- mcols(gi_coords)[, id_column]
        gi_coords$GRP_id <- ifelse(is.na(gi_coords$GRP_id), 0, gi_coords$GRP_id)
    } else {
        stop('id column: "', id_column, '" was not found in the gi object metadata')
    }

    # Use known seqlen for assembly hg38 or hg19, or require it to be provided
    if (genome == "hg38") {
        seql <- seqlengths(set_hg38_seqlengths(regions(gi_coords)))
    } else {
        if ((genome == "hg19")) {
            seql <- seqlengths(set_hg19_seqlengths(regions(gi_coords)))
        } else {
            # take from the regions
            seql <- seqlengths(regions(gi_coords))
        }
    }

    if (any(is.na(seql))) {
        stop("To output in sam file, seqlens are required. Please assign seqlengths
        to the acnhor regions of input gi object. Human hg38 and hg19 can be
        provided as 'genome' argument ")
    }

    temp_gi_coords <- gi_coords
    mcols(temp_gi_coords)$seq_anchor1 <- "*"
    mcols(temp_gi_coords)$seq_anchor2 <- "*"

    mcols(temp_gi_coords)$width_anchor1 <- width(temp_gi_coords)$first
    mcols(temp_gi_coords)$width_anchor2 <- width(temp_gi_coords)$second
    mcols(temp_gi_coords)$artQUAL1 <- "*"
    mcols(temp_gi_coords)$artQUAL2 <- "*"
    # create art CIGAR of perfect matches
    mcols(temp_gi_coords)$artCIGAR1 <- paste(mcols(temp_gi_coords)$width_anchor1, "M", sep = "")
    mcols(temp_gi_coords)$artCIGAR2 <- paste(mcols(temp_gi_coords)$width_anchor2, "M", sep = "")

    # there are 2 types of joints.
    # 1. chimeric junctions, relatively close and on SAME chromosome
    # 2. interchromosomal
    # get separated by cis and trans

    grange1 <- temp_gi_coords@regions[temp_gi_coords@anchor1]
    grange2 <- temp_gi_coords@regions[temp_gi_coords@anchor2]
    partner1_dist <- distance(grange1, grange2)
    mcols(temp_gi_coords)$partner1_dist <- partner1_dist

    message("Processing single-line sam records")
    cis_inter <- temp_gi_coords[!is.na(partner1_dist)]
    cis_inter_large_dist <- cis_inter[mcols(cis_inter)$partner1_dist >= distance_chim_junction]
    cis_inter <- cis_inter[mcols(cis_inter)$partner1_dist < distance_chim_junction]
    if (length(cis_inter) > 0) {
        cis_inter_tib <- as_tibble(cis_inter)
        # important to have same gene to have them correct
        mcols(cis_inter)$artCIGAR <- paste(mcols(cis_inter)$width_anchor1, "M", mcols(cis_inter)$partner1_dist, "N", mcols(cis_inter)$width_anchor2, "M", sep = "")
        mcols(cis_inter)$artSEQ <- "*"
        mcols(cis_inter)$artQUAL <- "*"
        mcols(cis_inter)$col2 <- 0
        mcols(cis_inter)$col5 <- 0
        mcols(cis_inter)$col7 <- "*"
        mcols(cis_inter)$col8 <- 0
        mcols(cis_inter)$col9 <- 0
        mcols(cis_inter)$col11 <- "*"
        mcols(cis_inter)$tag <- paste(paste("RG:Z:", cis_inter$GRP_id, sep = ""),
            paste("DG:Z:", cis_inter$GRP_id, sep = ""),
            paste("MR:Z:", paste(cis_inter_tib$seqnames2, cis_inter_tib$start2, cis_inter_tib$end2, cis_inter_tib$strand2, sep = ","), sep = ""),
            sep = "\t"
        )

        cis_inter_tibble <- as_tibble(cis_inter) %>%
            dplyr::select(
                rname, col2, seqnames1, start1,
                col5, artCIGAR, col7, col8, col9, artSEQ, col11, tag
            )
    } else {
        cis_inter_tibble <- tibble()
    }


    ###
    message("Processing double-line sam records")
    trans_inter <- temp_gi_coords[is.na(partner1_dist)]
    trans_inter <- c(trans_inter, cis_inter_large_dist)

    trans_inter_tib <- as_tibble(trans_inter)

    mcols(trans_inter)$col1read1 <- paste(mcols(trans_inter)[, "rname"], "", sep = "")
    mcols(trans_inter)$col1read2 <- paste(mcols(trans_inter)[, "rname"], "", sep = "")
    mcols(trans_inter)$col2 <- 0
    mcols(trans_inter)$col5 <- 0
    mcols(trans_inter)$col7 <- "*"
    mcols(trans_inter)$col8 <- 0
    mcols(trans_inter)$col9 <- 0
    mcols(trans_inter)$col11 <- "*"
    mcols(trans_inter)$tag <- paste(paste("RG:Z:", trans_inter$GRP_id, sep = ""),
        paste("DG:Z:", trans_inter$GRP_id, sep = ""),
        paste("MR:Z:", paste(trans_inter_tib$seqnames2,
            trans_inter_tib$start2,
            trans_inter_tib$end2,
            trans_inter_tib$strand2,
            sep = ","
        ), sep = ""),
        sep = "\t"
    )

    # as_tibble(trans_inter) %>% dplyr::select(duplex_id) %>% dplyr::summarise()
    # first I should have first pair of gi_coords
    tibble_read1 <- as_tibble(trans_inter) %>% dplyr::select(
        col1read1, col2,
        seqnames1, start1, col5,
        artCIGAR1, col7, col8,
        col9, seq_anchor1,
        col11, tag
    )
    tibble_read2 <- as_tibble(trans_inter) %>% dplyr::select(
        col1read2, col2,
        seqnames2, start2, col5,
        artCIGAR2, col7, col8,
        col9, seq_anchor2,
        col11, tag
    )

    ################################################################################
    # starting creating SAM fata.frame
    # header
    header_chr_sizes <- c("@HD\tVN:1.4\tSO:coordinate")
    # sizes from genome fasta
    # @SQ  SN:chr1  LN:248956422
    header_chr_sizes <- c(header_chr_sizes, paste("@SQ\tSN:", names(seql), "\tLN:", unname(seql), sep = ""))
    header_chr_sizes <- c(header_chr_sizes, "@RG\tID:sample1\tSM:sample1")
    header_chr_sizes <- c(header_chr_sizes, "@PG\tID:example_program\tPN:example_program\tVN:1.0")

    tibble_read1$seq_anchor1 <- unlist(tibble_read1$seq_anchor1)
    tibble_read2$seq_anchor2 <- unlist(tibble_read2$seq_anchor2)

    writeLines(text = header_chr_sizes, con = file_out)
    write.table(cis_inter_tibble, file_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
    write.table(tibble_read1, file_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
    write.table(tibble_read2, file_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)

    message("SAM is written into: ", file_out)
}
##########################################################################################
#
# library(BSgenome)
# gi = RNADuplexSampleClustReads
#
#
# #chr22:43,166,652-43,167,521
# gr_region = GRanges(IRanges(start=17504924,end=17524195),seqnames='chr22',strand = '+')
# sam_test = subsetByOverlaps(gi,gr_region)
# sam_test[sam_test$]
#
# sam_test = gi
# sam_cis = sam_both[sam_both@elementMetadata$cis==1]
#
# message("Load GENOME")
# #reference GENOME is required
# reference = "/data/meyer/egor/duplex_workflow/chr22.fa"
# genome <- readDNAStringSet(reference)
#
# gi = res$gi_reads
# writeGIToSAMfile(gi_coords=gi,genome,file_out = '../testsam.sam',distance_chim_junction = 1e5)
#
# output_sam_cis = paste("SamCis/",index,".sam",sep="")
# distance_chim_junction = 1000
#
# #classifyTwoArmChimeras(sam_test,min_junction_len = 10,junctions_gr = SampleSpliceJncGR,max_sj_shift = 5)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(SGSeq)
# txdb_v44 <- makeTxDbFromGFF("/fast/AG_Meyer/esemenc/SequenceData/hg38_transcriptome/gencode.v44.annotation.gtf")
# saveDb(txdb_v44,"../txdb_v44.txdb")
# tx_gtf<-AnnotationDbi::loadDb("../txdb_v44.txdb")
# txf <- convertToTxFeatures(tx_gtf)
# SampleSpliceJncGR<-GRanges(txf[txf@type=='J'])
#
# refjnc = subsetByOverlaps(SampleSpliceJncGR,gr_region+1e4)
# refjnc2 = subsetByOverlaps(txf[txf@type=='J'],gr_region+1e4)
# chimjnc = get_chimeric_junctions_onestrand(sam_test)
#
# chimjnc
# refjnc
#
# findOverlaps(chimjnc,jnc,maxgap = 5,type='equal')
#
# subsetByOverlaps(SampleSpliceJncGR,get_chimeric_junctions_onestrand(sam_test),type ='equal',maxgap=10)


# library(AnnotationHub)
# library(Biostrings)
#
# # Create an AnnotationHub object
# ah <- AnnotationHub()
#
# # Query AnnotationHub for the hg38 chromosome 22 sequence
# query_result <- query(ah, c("Homo sapiens",'fa','2bit'))
#
# ah2 <- query(ah, c("fasta", "homo sapiens", "Ensembl",'primary_assembly'))
# ah2[[ah2[1]$ah_id]]
# s = readDNAStringSet('../../chr22.fa')
# rtracklayer::export.2bit(con = '../../chr22.2bit', test_2bit_out)
#
#
# # Find the specific resource for the chromosome 22 FASTA file
# # Here we assume the first hit is the desired one for simplicity, but you might want to check `query_result` for the exact resource
# chr22_resource <- query_result[[1]]
#
# # Download the resource
# chr22_fasta <- chr22_resource$resource
# ah[query_result[1]$ah_id]
# # Read the FASTA file into a DNAStringSet object
# chr22_dna <- readDNAStringSet(query_result[1]$ah_id)
#
# # Return the DNAStringSet object
# return(chr22_dna)
