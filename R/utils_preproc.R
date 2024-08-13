#' Check the column names and types in read dataframe
#' @description
#'
#' Function to check the correct column names and types in dataframe input.
#' Tries to guess the column names, if colnames are not provided,
#' but the types are correct
#'
#' @details
#' Expected column names for bedpe file
#' \code{ c("chromA","startA",'endA',"chromB",
#'                     'startB','endB','readname','flag','strandA','strandB')}
#'  Expected colnames for STAR Chimeric junction input
#'  For the 'old' chimeric detection scheme
#' \code{ c("chr_donorA","brkpt_donorA",
#' "strand_donorA","chr_acceptorB",
#' "brkpt_acceptorB","strand_acceptorB",
#' "junction_type","repeat_left_lenA",
#' "repeat_right_lenB","read_name",
#' "start_alnA","cigar_alnA",
#' "start_alnB","cigar_alnB",
#' "num_chim_aln","max_poss_aln_score",
#' "non_chim_aln_score","this_chim_aln_score",
#' "bestall_chim_aln_score") }
#' For the 'new' chimeric detection scheme
#' \code{c("chr_donorA","brkpt_donorA",
#'                   "strand_donorA","chr_acceptorB",
#'                   "brkpt_acceptorB","strand_acceptorB",
#'                   "junction_type","repeat_left_lenA",
#'                   "repeat_right_lenB","read_name",
#'                   "start_alnA","cigar_alnA",
#'                   "start_alnB","cigar_alnB")}
#'
#' @keywords internal
#' @param df input
#' @param table_type one in c("STAR","bedpe")
#'
#' @return dataframe with the properly formatted columns
col_check_rename <- function(df, table_type = "STAR") {
    ncol <- dim(df)[2]

    # load without storing globally
    STAR_COLNAMES19 <- get_colnames_and_types_for_input("STAR_COLNAMES19")
    STAR_COLNAMES14 <- get_colnames_and_types_for_input("STAR_COLNAMES14")
    STAR_COLTYPES19 <- get_colnames_and_types_for_input("STAR_COLTYPES19")
    STAR_COLTYPES14 <- get_colnames_and_types_for_input("STAR_COLTYPES14")
    BEDPE_COLNAMES <- get_colnames_and_types_for_input("BEDPE_COLNAMES")
    BEDPE_COLTYPES <- get_colnames_and_types_for_input("BEDPE_COLTYPES")


    if (table_type != "STAR") {
        expected_names <- BEDPE_COLNAMES
        expected_types <- BEDPE_COLTYPES
        n_col_check <- 10
        message("Detected input type generic .bedpe dataframe , expecting 10+ columns \n")
    }
    if (table_type == "STAR") {
        if (all(STAR_COLNAMES19 %in% colnames(df))) {
            expected_names <- STAR_COLNAMES19
            expected_types <- STAR_COLTYPES19
            n_col_check <- 19
            message("Detected input type Chimeric.out.Junction with multimaps, expecting first ", n_col_check, "columns to follow STAR naming\n")
        } else {
            if (all(STAR_COLNAMES14 %in% colnames(df))) {
                expected_names <- STAR_COLNAMES14
                expected_types <- STAR_COLTYPES14
                message("Detected input type Chimeric.out.Junction without multimaps, expecting first ", n_col_check, "columns  to follow STAR naming \n")
                n_col_check <- 14
            } else {
                expected_names <- STAR_COLNAMES14
                expected_types <- STAR_COLTYPES14
                n_col_check <- 14
                message("Specified format is STAR Chimeric.out.Junction, but columns are missing or unnamed \n")
                message("Will guess and validate first 14 columns. Consider transformig input to more general .bedpe format or settignt he column names ")
            }
        }
    }

    df <- as_tibble(df)
    testt <- data.frame(head(df))

    actual_names <- names(testt)
    actual_types <- vapply(
        testt %>% dplyr::select(all_of(seq_len(n_col_check))),
        function(col) class(col),""
    )

    # Check if there are extra columns
    extra_columns <- setdiff(actual_names, expected_names)

    # check types
    if (all(actual_types == expected_types)) {
        message("Verified the correct column types and order \n")
        # if(length(extra_columns) > 0) {
        #   message("Extra columns are found:", paste(extra_columns, collapse = ", "), "\n")
        # } else {
        #   message("No extra columns were found.\n")
        # }
        # check names
        if (!all(expected_names %in% colnames(testt))) {
            message("Guessing columns")
            df <- df %>%
                rename_with(~expected_names, seq_len(n_col_check))
            status <- data.frame(Before = actual_names[seq_len(n_col_check)], After = names(df)[seq_len(n_col_check)])
            message(paste(status$Before, "\t-> ", status$After, collapse = "\n", sep = ""))
            message("\n ")
        } else {
            message("Column names are correct \n")
        }
    } else { # TYPES INCONSISTENCY

        result_guess <- tryCatch(
            {
                df_tmp <- head(df, 5)
                df_tmp[seq_len(n_col_check)] <- lapply(seq_len(n_col_check), function(i) as(df_tmp[[i]], expected_types[i]))
                TRUE # Return TRUE if successful
            },
            error = function(e) {
                FALSE # Return FALSE if conversion of type is not possible
            }
        )
        # if types can be formatted, do that and check names
        if (result_guess) {
            df[seq_len(n_col_check)] <- lapply(seq_len(n_col_check), function(i) as(df[[i]], expected_types[i]))
            df <- df %>%
                rename_with(~expected_names, all_of(actual_names[seq_len(n_col_check)]))
            status <- data.frame(Before = actual_names[seq_len(n_col_check)], After = names(df)[seq_len(n_col_check)])
            message(paste(status$Before, "\t-> ", status$After, collapse = "\n", sep = ""))
            msgtext <- paste0(colnames(df), sep = " ")
            message("\nGuessed column names and converted types \n", msgtext)


        } else {
            # if something was wrong with types guessing
            stats <- ifelse(actual_types[seq_len(n_col_check)] == expected_types, "ok", "<-problem")
            expected <- paste(seq_len(n_col_check),
                expected_names, ":\t",
                expected_types, "\t|| ",
                actual_names[seq_len(n_col_check)], ":\t", actual_types[seq_len(n_col_check)], "\t:",
                stats,
                collapse = "\n"
            )
            message("Problem with the input format found")
            message("Expected name :\t type \t || Found name : type\n", expected, sep = "")
            stop("\nColumn names, types, or order do not match the expected criteria.")
        }
    }
    return(df)
}

#' Processing of of the STAR SE Chimeric.junction.out
#'
#' Calculates alignment coordinates and returns reads with categories
#'
#' #' \describe{
#'   \item{multimap}{ multi-mapped read}
#'   \item{multigap}{ more than one junction (more than two 'N' in CIGAR string)}
#'   \item{bad junction}{ Artifacts. I.e alignments for both arms are continious, but with 'backward' chimeric junction was wrongly put}
#' }
#' @keywords interal
#' @param dt Chimeric.out.junction with the correct column names
#' @param keep_all_columns - TRUE or FALSE. Keep CIGAR strings and junction coordinate columns
#' @seealso [col_check_rename()]
#' @return tibble with annotated reads
preproc_chim_junction_out_se <- function(dt, keep_all_columns = FALSE) {
    message("processing raw Chimeric.junction.out...")

    if (!("this_chim_aln_score" %in% colnames(dt))) {
        dt$this_chim_aln_score <- 1
        dt$num_chim_aln <- 1
    }
    if (!keep_all_columns) {
        dt <- dt %>% dplyr::select(c(
            read_name,
            chr_donorA, chr_acceptorB,
            start_alnA, start_alnB,
            strand_donorA, strand_acceptorB,
            cigar_alnA, cigar_alnB,
            brkpt_donorA, brkpt_acceptorB, this_chim_aln_score,
            num_chim_aln
        ))
    }
    # rename, count gaps
    dt <- dt %>%
        dplyr::rename(
            readname = read_name,
            chromA = chr_donorA,
            strandA = strand_donorA,
            startA = start_alnA,
            chromB = chr_acceptorB,
            strandB = strand_acceptorB,
            startB = start_alnB,
            score = this_chim_aln_score
        ) %>%
        mutate(
            ngapsA = str_count(cigar_alnA, "N"),
            ngapsB = str_count(cigar_alnB, "N")
        ) %>%
        mutate(
            lengapsA = get_char_count_cigar(cigar_alnA, "N"),
            lengapsB = get_char_count_cigar(cigar_alnB, "N")
        )

    message("getting alignment coordinates ")
    dt <- dt %>%
        mutate(
            endA = startA + GenomicAlignments::cigarWidthAlongReferenceSpace(cigar_alnA),
            endB = startB + GenomicAlignments::cigarWidthAlongReferenceSpace(cigar_alnB),
            bad_junction = ifelse(chromA == chromB &
                strandA == strandB &
                abs(brkpt_donorA - brkpt_acceptorB) <= 2, 1, 0),
            bad_junction = ifelse(
                chromA == chromB &
                    strandA == strandB &
                    (abs(startA - endB) <= 3 | abs(endA - startB) <= 3), 1, bad_junction
            ),
            multigap = as.integer((lengapsA >= 2) | (lengapsB >= 2)),
            multimap = as.integer(num_chim_aln > 1)
        ) %>%
        relocate(endB, .after = startB) %>%
        relocate(endA, .after = startA)

    if (!keep_all_columns) {
        dt <- dt %>% dplyr::select(-c(
            brkpt_donorA, brkpt_acceptorB,
            cigar_alnA, cigar_alnB, num_chim_aln
        ))
    }
    return(dt)
}

#' Processing of of the STAR PE Chimeric.junction.out
#'
#' Calculates alignment coordinates and returns reads with categories
#'
#' #' \describe{
#'   \item{multimap}{ multi-mapped read}
#'   \item{multigap}{ more than one junction (more than two 'N' in CIGAR string)}
#'   \item{bad junction}{ Artifacts. I.e alignments for both arms are continious, but with 'backward' chimeric junction was wrongly put}
#' }
#' @keywords internal
#' @param dt Chimeric.out.junction with the correct column names
#' @param keep_all_columns - TRUE or FALSE. Keep CIGAR strings and junction coordinate columns
#'
#' @return tibble with annotated reads
preproc_chim_junction_out_pe <- function(dt, keep_all_columns = FALSE) {
    message("processing raw Chimeric.junction.out... PE")
    bad_junction_tr = 6
    if (!keep_all_columns) {
        dt <- dt %>% dplyr::select(c(
            read_name,
            chr_donorA, chr_acceptorB,
            start_alnA, start_alnB,
            strand_donorA, strand_acceptorB,
            cigar_alnA, cigar_alnB,
            brkpt_donorA, brkpt_acceptorB, this_chim_aln_score,
            num_chim_aln
        ))
    }

    # rename, count gaps
    message("getting alignment coordinates ")
    dt <- dt %>%
        dplyr::rename(
            readname = read_name,
            chromA = chr_donorA,
            strandA = strand_donorA,
            startA = start_alnA,
            chromB = chr_acceptorB,
            strandB = strand_acceptorB,
            startB = start_alnB,
            score = this_chim_aln_score
        ) %>%
        mutate(
            ngapsA = str_count(cigar_alnA, "N"),
            ngapsB = str_count(cigar_alnB, "N")
        ) %>%
        add_column(
            lengapsA = 0,
            lengapsB = 0
        ) %>%
        mutate(
            cigar_str = cigar_alnA,
            has_p = str_count(cigar_str, "p"),
            minus_p = str_count(cigar_str, "-"),
            sign_p = if_else(minus_p == 1, -1, 1 * has_p),
            n_p = get_char_count_cigar(cigar_str, "p"),
            n_s = get_char_count_cigar(cigar_str, "S"),
            n_m = get_char_count_cigar(cigar_str, "M"),
            n_i = get_char_count_cigar(cigar_str, "I"),
            n_n = get_char_count_cigar(cigar_str, "N"),
        ) %>%
        mutate(
            lengapsA = n_n,
            endA_along = startA + cigarWidthAlongReferenceSpace(cigar_str),
            lA =  n_m + n_p * sign_p + n_n,
            endA = startA + lA - 1
        ) %>%
        mutate(
            cigar_str = cigar_alnB,
            has_p = str_count(cigar_str, "p"),
            minus_p = str_count(cigar_str, "-"),
            sign_p = if_else(minus_p == 1, -1, 1 * has_p),
            n_p = get_char_count_cigar(cigar_str, "p"),
            n_s = get_char_count_cigar(cigar_str, "S"),
            n_m = get_char_count_cigar(cigar_str, "M"),
            n_i = get_char_count_cigar(cigar_str, "I"),
            n_n = get_char_count_cigar(cigar_str, "N")
        ) %>%
        mutate(
            lengapsB = n_n,
            lB = n_m + n_p * sign_p + n_n,
            endB_along = startB + cigarWidthAlongReferenceSpace(cigar_str),
            endB = startB + lB - 1
        ) %>%
        mutate(
            bad_junction = ifelse(chromA == chromB &
                strandA == strandB &
                (abs(startA - endB) <= 3) |
                (abs(endA - startB) <= 3), 1, 0),
            multigap = as.integer((lengapsA >= 2) | (lengapsB >= 2)),
            multimap = as.integer(num_chim_aln > 1)
        ) %>%
        relocate(endB, .after = startB) %>%
        relocate(endA, .after = startA) %>%
        dplyr::select(-c(
            cigar_str, has_p, minus_p, sign_p,
            n_p, n_s, n_m, n_i, n_n,
            ngapsA, ngapsB, lengapsA, lengapsB,lA, lB
        ))

    if (!keep_all_columns) {
        dt <- dt %>% dplyr::select(-c(
            brkpt_donorA, brkpt_acceptorB,
            cigar_alnA, cigar_alnB, num_chim_aln
        ))
    }

    return(dt)
}

#' Preprocess .bedpe input
#'
#' Searches for the multi-mapped and bad reads (overlapping arms)
#' Adds 'multimap', 'bad_junction' columns filled with 0 or 1 and 'multigap' = 0
#' for consistency with other pre-processing methods
#'
#' @keywords internal
#' @param dt dataframe with reads aligned to strictly two loci
#' @param keep_all_columns keeep columns apart form the required from .bedpe format
#'
#' @return pre-processed dataframe
preproc_generic <- function(dt, keep_all_columns = TRUE) {
    message("processing generic bedpe input...")
    BEDPE_COLNAMES <- get_colnames_and_types_for_input("BEDPE_COLNAMES")
    if (!keep_all_columns) {
        dt <- dt %>% dplyr::select(all_of(BEDPE_COLNAMES))
    }
    dt <- dt %>% mutate(
        bad_junction = 0,
        multimap = 0,
        multigap = 0
    )
    dt <- dt %>%
        mutate(
            bad_junction = ifelse(chromA == chromB & strandA == strandB &
                (abs(startA - endB) <= 3 | abs(endA - startB) <= 3), 1, bad_junction)
        ) %>%
        group_by(readname) %>%
        mutate(
            nr = n(),
            multimap = ifelse(nr != 1, 1, multimap)
        ) %>%
        ungroup()
    dt$nr <- NULL
    if (!("score" %in% colnames(dt))) {
        message("no alignment score column found. Will create one for consistency")
        dt$score <- 1
    }


    return(dt)
}

#' Preprocess GInteractions input
#'
#' Searches for the multi-mapped reads (overlapping arms)
#' Adds 'multimap', 'bad_junction' columns filled with 0/1 and 'multigap' = 0
#' for consistency with other pre-processing methods.
#'
#' @keywords internal
#' @param gi_raw `GInteractions` with inpit RNA interactions
#' @param keep_all_columns keep columns apart from those required by .bedpe format
#'
#' @return pre-processed dataframe
preproc_generic_gi <- function(gi_raw, keep_all_columns = TRUE) {
    message("processing generic GInteractions input...")
    BEDPE_COLNAMES <- get_colnames_and_types_for_input("BEDPE_COLNAMES")
    cnames <- colnames(mcols(gi_raw))
    dt <- makeDfFromGi(gi_raw)
    if (!("readname" %in% cnames)) {
        message("no radname column found. Will create one for consistency")
        dt$readname <- seq_along(gi_raw)
    }
    if (!("score" %in% cnames)) {
        message("no alignment score column found. Will create one for consistency")
        dt$score <- 1
    }
    if (!keep_all_columns) {
        dt <- dt %>% dplyr::select(any_of(BEDPE_COLNAMES), score)
    }
    dt <- dt %>% mutate(
        bad_junction = 0,
        multimap = 0,
        multigap = 0
    )
    dt <- dt %>%
        mutate(
            bad_junction = ifelse(chromA == chromB & strandA == strandB &
                (abs(startA - endB) <= 3 | abs(endA - startB) <= 3), 1, bad_junction)
        ) %>%
        group_by(readname) %>%
        mutate(
            nr = n(),
            multimap = ifelse(nr != 1, 1, multimap)
        ) %>%
        ungroup()
    dt$nr <- NULL

    return(dt)
}


#' Run pre-processing of chimeric reads input
#' @description
#' Imports dataframe with reads (.bedpe or Chimeric.out.junction ) or `GInteractions`
#' object. Adds necessary metadata
#' If the input is
#' Checks column names or tries to quess them if not provided.
#' For STAR input, calculates length of the alignment and annotates alignments types
#'
#' @details
#' If not existed, adds fields required for the downstream steps:
#' 'readname', 'map_type', 'score', 'n_reads'.
#' 'map_type' field determines the type of the chimeric read:
#' \describe{
#'   \item{multimap}{ multi-mapped read}
#'   \item{multigap}{ more than one junction (more than two 'N' in CIGAR string)}
#'   \item{bad junction}{ Artifacts or possibly unaccounted types.
#'   I.e alignments for both arms are continuous, but with 'backward' chimeric
#'   junction was wrongly introduced in the mapping}
#' }
#' @param data Either dataframe-like object: Chimeric.out.junction from STAR or
#' .bepde - formatted or `GInteractions` object from \pkg{InteractionSet} package
#' @param keep_metadata - TRUE or FALSE. do not extra fields like
#'  CIGAR strings and junction coordinates
#' @param table_type in \code{c("STAR","bedpe")} for Chimeric.out.Junction or generic input
#' @param library_type \code{c("SE","PE")} for pair- or single- end input
#' @param return_gi if the return object should be `GInteractions`
#' @param  min_arm_len minimum allowed length of the alignment arm.
#' Read will be dropped if either arm is shorter
#' @return tibble with new metadata fields OR GInteractions if `return_gi` is
#' set to TRUE
#' @export
#' @examples
#' # load data
#' data(RNADuplexesSampleData)
#' # with bedpe input
#' preproc_reads <- runDuplexDiscoPreproc(RNADuplexesRawBed, table_type = "bedpe")
#' # with STAR input
#' preproc_reads_star <- runDuplexDiscoPreproc(RNADuplexesRawChimSTAR,
#'     table_type = "STAR",
#'     keep_metadata = FALSE
#' )
runDuplexDiscoPreproc <- function(data, table_type,
    library_type = "SE",
    keep_metadata = TRUE,
    return_gi = FALSE,
    min_arm_len = 15) {
    # Start with determining the input data type
    if (is(data, "GInteractions")) {
        # GI input
        message("Detected input type is GInteractions")
        df <- preproc_generic_gi(gi_raw = data, keep_all_columns = keep_metadata)
    } else {
        # DF input
        if (is(data, "data.frame")) {
            message("Input type is table")
            if (!table_type %in% c("STAR", "bedpe")) {
                stop("Input type of table is not specified or not accepted. Accepted values are: bedpe, STAR ")
            }
            # Table input preproc
            if (!(table_type %in% c("STAR", "bedpe"))) {
                message("Wrong type of the table provided. Use 'STAR' for
              Chimeric.out.Junction or \n 'bedpe' for generic bedpe format \n")
            }
            df <- data
            df <- col_check_rename(df, table_type)

            if (table_type == "STAR" && (library_type %in% c("SE", "PE"))) {
                if (library_type == "SE") {
                    df <- preproc_chim_junction_out_se(df, keep_all_columns = keep_metadata)
                }
                if (library_type == "PE") {
                    df <- preproc_chim_junction_out_pe(df, keep_all_columns = keep_metadata)
                }
            } else {
                df <- preproc_generic(df, keep_all_columns = keep_metadata)
            }
        } else {
            stop("cannot determine input data type. Provided : ", class(data)[1])
        }
    }
    
    # Add alignment length 
    df <- df %>% mutate(aln_lenA = endA - startA,
                         aln_lenB = endB - startB)
  
    # Add number of reads per record (should be 1, but check in case of collapsed input)
    if (!("n_reads" %in% colnames(df))) {
        df <- df %>% mutate(
            n_reads = 1,
            read_id = seq_len(nrow(df))
        )
    } else {
        df <- df %>% add_column(read_id = seq_len(nrow(df)))
    }
    df = df %>%  mutate(
        len_too_short = as.integer(!(aln_lenA > min_arm_len &
                                     aln_lenB > min_arm_len)))
    
    df <- df %>% mutate(map_type = case_when(
        multigap == 0 & multimap == 1 & bad_junction == 0 ~ "multi_map",
        multigap == 1 & multimap == 0 & bad_junction == 0 ~ "multi_split",
        multigap == 1 & multimap == 1 & bad_junction == 0 ~ "multi_split&map",
        bad_junction == 1 ~ "not_chim",
        len_too_short == 1 ~ "too_short",
        multigap == 0 & multimap == 0 & bad_junction == 0 ~ "2arm",
        .default = "notype"
    ))

    df <- df %>%
        dplyr::select(-c(multigap, multimap, bad_junction)) %>%
        dplyr::relocate(c(chromA, startA, endA, strandA, chromB, startB, endB, strandB, readname, map_type, score), .before = everything())
    if (!keep_metadata) {
        df <- df %>% dplyr::select(c(chromA, startA, endA, strandA, chromB, startB, endB, strandB, readname, map_type, score, n_reads))
    }
    if (return_gi) {
        df <- makeGiFromDf(df)
    }
    return(df)
}


.filterAlignmentLength <- function(dt,minlen){
  
  dt = dt %>% dplyr::mutate(len_not_ok=!(aln_lenA > minlen & aln_lenB > minlen))
  nf = sum(as.integer(dt$len_not_ok))
  nr = nrow(dt)
  dt  = dt %>% dplyr::filter(!len_not_ok)
  dt$len_not_ok = NULL
  message("Filtered out : ",nf,
          " too short alignments out of ", nr," : ", round(nf/nr*100,2) , ' %')
  return(dt)
}

# accepts single gap df
# integrate after preproc, and after SJ ? 
#' Extract regions around chimeric junction 
#'  
#' Trim alignements to contain only 'extract len' nucleotides adajcent
#' to the chimeric junction
#'  
#' @details
#' In case of the long alignemtns, it may be necessary trim chimeric alignments
#' to identify RNA duplex. If 'extract_len' is longer than the read 
#' alignemnt length, then no trimmin is performed
#'   
#' @param dt table with the 
#' @param extract_len 
#'
#' @return dataframe with the trimmed alignments 
#' @export
#' @examples
#' data("RNADuplexesSampleData")
#' dt_preproc = runDuplexDiscoPreproc(RNADuplexesRawChimSTAR,
#' table_type = 'STAR',library_type = 'SE')
#' trimAroundJunction(dt_preproc,40)
#' 
trimAroundJunction <- function(dt,
                               extract_len = 30){
  cnames = c('brkpt_donorA','brkpt_acceptorB')
  dift = 2 
  extract_len = extract_len -1
  
  if (!all(cnames %in% colnames(dt))){
    message("cannot extract region around junction
            because junction fields are not present")
    return(dt)
  }
  
  # Types of junction arrangement (point of ligation)
  # t1: normal ==* *==
  # t2: backwards *== ==*
  # t3: A normal B backwards  ==* ==*
  # t2: B normal A backwards *== *==
  if (all(c("jA","jB") %in% colnames(dt))){
    #pass
  }else{
    if (all(cnames %in% colnames(dt))){
      dt <- dt %>% rename(jA = brkpt_donorA,jB = brkpt_acceptorB )
    }else{
      message("Junction fields not found: Possible filed names are jA and JB or 
             brkpt_donorA and brkpt_acceptorB ")
    }
  }
  
  dt1 = dt  %>%
    mutate(t1 = if_else( (abs(endA - jA) <= dift ) & (abs(startB - jB) <= dift),1,0 ),
           t2 = if_else( (abs(startA - jA)<= dift) & (abs(endB - jB) <= dift) ,1,0),
           t3 = if_else( (abs(endA - jA)<= dift) & (abs(endB - jB) <= dift) ,1,0),
           t4 = if_else( (abs(startA - jA)<= dift) & (abs(startB - jB) <= dift) ,1,0),
           any1 = t1+t2+t3+t4) 
  
  
  if (!any(dt1$any1>1) & !any(dt1$any1==0)){
    message('type detected')
  }else{
    message('some problem occured in detectin jucntion type.
            Check alignment lengths')
  }
  
  dt2 = dt1 %>%
    tidyr::pivot_longer(
      cols = all_of(c('t1','t2','t3','t4')),
      names_to = "junctype",
      values_to = "value"
    ) %>%
    filter(value == 1) %>%
    select(-value)
  

   dt_modified <- dt2 %>%
    mutate(
      aln_lenA = endA - startA,
      aln_lenB = endB - startB,
      extrA = pmin(aln_lenA, extract_len),
      extrB = pmin(aln_lenB, extract_len)) %>%
    mutate(
      startA = case_when(
        junctype == 't1' ~ endA - extrA,
        junctype == 't2' ~ startA,
        junctype == 't3' ~ endA - extrA,
        junctype == 't4' ~ startA
      ),
      endA = case_when(
        junctype == 't1' ~ endA,
        junctype == 't2' ~ startA + extrA,
        junctype == 't3' ~ endA,
        junctype == 't4' ~ startA + extrA
      ),
      startB = case_when(
        junctype == 't1' ~ startB,
        junctype == 't2' ~ endB - extrB,
        junctype == 't3' ~ endB - extrB,
        junctype == 't4' ~ startB
      ),
      endB = case_when(
        junctype == 't1' ~ startB + extrB,
        junctype == 't2' ~ endB,
        junctype == 't3' ~ endB,
        junctype == 't4' ~ startB + extrB
      )
    )
  #dt_modified$junctype = NULL
  return(dt_modified)
}

