#' Pairs parsing functions
#' 
#' @param file pairs file. Default formatting is `<readname>\t<chr1>\t<start1>\t<chr2>\t<start2>`. 
#' @param chr1.field,start1.field,chr2.field,start2.field,strand1.field,strand2.field,frag1.field,frag2.field Index of the column in which each field is contained in the pairs file.
#' @param nThread Number of CPUs to use to import the `pairs` file in R
#' @param nrows Number of pairs to import
#' @return a GInteractions object
#'
#' @importFrom vroom vroom
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom IRanges IRanges
#' @rdname parse-pairs
#' @keywords internal

.pairs2gi <- function(
    file, 
    chr1.field = NULL, 
    start1.field = NULL, 
    chr2.field = NULL, 
    start2.field = NULL, 
    strand1.field = NULL, 
    strand2.field = NULL, 
    frag1.field = NULL, 
    frag2.field = NULL, 
    nThread = 1, 
    nrows = Inf  
) {

    # -- Guess pairs format
    skip <- 0
    pair_type.field <- NULL
    if (
        is.null(chr1.field) & 
        is.null(start1.field) & 
        is.null(chr2.field) & 
        is.null(start2.field) & 
        is.null(strand1.field) & 
        is.null(strand2.field) & 
        is.null(frag1.field) & 
        is.null(frag2.field)
    ) {
        .fmt <- .guess_pairs_format(file)
    } 
    else {
        .fmt <- 'prespecified'
    }
    if (.fmt == '4dn') {
        fileComments <- read.delim(
            file,
            sep = '&',
            nrows = 1000,
            check.names = FALSE, 
            header = FALSE 
        ) 
        fileComments <- fileComments[grepl('^#[^#]', fileComments[,1]), ]
        if (any(grepl("^#columns:", fileComments))) {
            colIDs <- fileComments[which(grepl("^#columns:", fileComments))] |> 
                strsplit(" ") |> 
                unlist() |> 
                tail(-1)
            chr1.field <- which(colIDs == 'chr1')
            start1.field <- which(colIDs == 'pos1')
            chr2.field <- which(colIDs == 'chr2')
            start2.field <- which(colIDs == 'pos2')
            strand1.field <- which(colIDs == 'strand1')
            strand2.field <- which(colIDs == 'strand2')
            if ("pair_type" %in% colIDs) {
                pair_type.field <- which(colIDs == 'pair_type')
            } else {
                pair_type.field <- NULL
            }
            if ("frag1" %in% colIDs) frag1.field <- which(colIDs == 'frag1')
            if ("frag2" %in% colIDs) frag2.field <- which(colIDs == 'frag2')
        }
        else {
            chr1.field <- 2
            start1.field <- 3
            chr2.field <- 4
            start2.field <- 5
            strand1.field <- 6
            strand2.field <- 7
        }
    }
    if (.fmt == 'hicpro') {
        # <ID> <chr1> <pos1> <str1> <chr2> <pos2> <str2> <isize> <frag1> <frag2> <mapq1> <mapq2> [<allele-spe. info>]
        fileLines <- vroom::vroom(
            file,
            n_max = 1000,
            col_names = FALSE,
            comment = '#',
            progress = FALSE, 
            show_col_types = FALSE
        )
        chr1.field <- 2
        start1.field <- 3
        chr2.field <- 5
        start2.field <- 6
        strand1.field <- 4
        strand2.field <- 7
        if (ncol(fileLines) >= 10) {
            frag1.field <- 8
            frag2.field <- 9
        }
    }
    if (.fmt == 'hicstuff') {
        skip <- 1
        # <chr1> <pos1> <end1> <chr2> <start2> <end2> <frag1> <frag2> <kernel_id> <iteration> <score> <pvalue> <qvalue>
        fileLines <- vroom::vroom(
            file,
            skip = 1, 
            n_max = 1000,
            col_names = FALSE,
            comment = '#',
            progress = FALSE, 
            show_col_types = FALSE
        )
        chr1.field <- 1
        start1.field <- 2
        chr2.field <- 4
        start2.field <- 5
        strand1.field <- 0
        strand2.field <- 0
        frag1.field <- 0
        frag2.field <- 0
    }
    if (any(c(
        is.null(chr1.field), 
        is.null(start1.field), 
        is.null(strand1.field), 
        is.null(chr2.field), 
        is.null(start2.field), 
        is.null(strand2.field)
    ))) stop("Impossible to parse pairs file")

    # -- Import pairs as anchors
    sel1 <- (c(chr1.field, start1.field, strand1.field, frag1.field))
    sel1 <- sel1[sel1 != 0]
    sel2 <- (c(chr2.field, start2.field, strand2.field, frag2.field))
    sel2 <- sel2[sel2 != 0]
    sel <- c(sel1, sel2)
    anchors <- vroom::vroom(
        file,
        skip = skip,
        n_max = nrows,
        col_select = dplyr::all_of(sel),
        comment = '#',
        col_names = FALSE,
        show_col_types = FALSE,
        num_threads = nThread
    )
    anchors1 <- dplyr::select(anchors, seq_len(ncol(anchors)/2))
    anchors2 <- dplyr::select(anchors, ncol(anchors)/2 + seq_len(ncol(anchors)/2))
    anchor_one <- GenomicRanges::GRanges(
        anchors1[[1]],
        IRanges::IRanges(anchors1[[2]], width = 1)
    )
    anchor_two <- GenomicRanges::GRanges(
        anchors2[[1]],
        IRanges::IRanges(anchors2[[2]], width = 1)
    )

    # -- Add strand information for anchors
    if (strand1.field != 0) GenomicRanges::strand(anchor_one) <- anchors1[[3]]
    if (strand2.field != 0) GenomicRanges::strand(anchor_two) <- anchors2[[3]]

    # -- Parse anchors into GInteractions
    gi <- InteractionSet::GInteractions(anchor_one, anchor_two)

    # -- Add frag information
    if (!is.null(frag1.field) & !is.null(frag2.field)) {
        gi$frag1 <- anchors1[[ncol(anchors1)]]
        gi$frag2 <- anchors2[[ncol(anchors1)]]
    } 
    else {
        gi$frag1 <- NA
        gi$frag2 <- NA
    }

    # -- Add pair_type information
    if (!is.null(pair_type.field)) {
        pair_types <- vroom::vroom(
            file,
            skip = skip,
            n_max = nrows,
            col_select = dplyr::all_of(pair_type.field),
            comment = '#',
            col_names = FALSE,
            show_col_types = FALSE,
            num_threads = nThread
        )
        gi$pair_type <- pair_types[[1]]
    } 

    # -- Add loop scores for hicstuff files 
    if (.fmt == 'hicstuff') {
        dat <- vroom::vroom(
            file,
            skip = 0,
            n_max = nrows,
            col_select = c(9, 10, 11, 12, 13),
            comment = '#',
            col_names = TRUE,
            show_col_types = FALSE,
            num_threads = nThread
        )
        GenomicRanges::mcols(gi) <- cbind(GenomicRanges::mcols(gi), dat)
    }

    # -- Add pairdist
    gi$distance <- InteractionSet::pairdist(gi) 

    # -- Add seqinfo
    if (.fmt == '4dn') {
        chrs <- fileComments |> 
            grep('#chromsize: ', x = _, value = TRUE) |>
            gsub("#chromsize: ", "", x = _)
        si <- gsub(".* ", "", chrs)
        names(si) <- gsub(" .*", "", chrs)
        GenomeInfoDb::seqlevels(gi) <- names(si)
        GenomeInfoDb::seqinfo(gi) <- GenomeInfoDb::Seqinfo(
            names(si), 
            as.numeric(si)
        )
    }

    return(gi)
}

.guess_pairs_format <- function(file) {

    lines <- readLines(file, n = 1000)
    fileComments <- grep('^#', lines, value = TRUE)
    fileLines <- vroom::vroom(
        file,
        n_max = 1000,
        col_names = FALSE,
        comment = '#',
        progress = FALSE, 
        show_col_types = FALSE, 
        skip = 1
    )
    ncols <- ncol(fileLines)

    # -- Check if file is pairtools format (pairtools)
    # <ID> <chr1> <pos1> <chr2> <pos2> <str1> <str2> <pair_type>
    if (any(grepl('pairs format', fileComments))) return('4dn')
    if (tryCatch({
        ncol(fileLines) >= 7 & 
        is.character(fileLines[[1]]) & # <ID>
        {is.character(fileLines[[2]]) | is.numeric(fileLines[[2]])} & # <chr1> 
        is.numeric(fileLines[[3]]) & # <pos1> 
        {is.character(fileLines[[4]]) | is.numeric(fileLines[[4]])} & # <chr2> 
        is.numeric(fileLines[[5]]) & # <pos2> 
        is.character(fileLines[[6]]) & all(unique(fileLines[[6]]) %in% c('+', '-')) & # <str1> 
        is.character(fileLines[[7]]) & all(unique(fileLines[[7]]) %in% c('+', '-')) # <str2>
    }, error = function(e) FALSE)) return('4dn')

    # -- Check if file is .pairs format (4DN)
    # <ID> <chr1> <pos1> <chr2> <pos2> <str1> <str2> [<frag1> <frag2>]
    if (any(grepl('pairs format', fileComments))) return('4dn')
    if (tryCatch({
        ncol(fileLines) >= 7 & 
        is.character(fileLines[[1]]) & # <ID>
        {is.character(fileLines[[2]]) | is.numeric(fileLines[[2]])} & # <chr1> 
        is.numeric(fileLines[[3]]) & # <pos1> 
        {is.character(fileLines[[4]]) | is.numeric(fileLines[[4]])} & # <chr2> 
        is.numeric(fileLines[[5]]) & # <pos2> 
        is.character(fileLines[[6]]) & all(unique(fileLines[[6]]) %in% c('+', '-')) & # <str1> 
        is.character(fileLines[[7]]) & all(unique(fileLines[[7]]) %in% c('+', '-')) # <str2>
    }, error = function(e) FALSE)) return('4dn')

    # -- Check if file is .validPairs format (HiC-Pro)
    # <ID> <chr1> <pos1> <str1> <chr2> <pos2> <str2> <isize> <frag1> <frag2> <mapq1> <mapq2> [<allele-spe. info>]
    if (tryCatch({
        ncol(fileLines) >= 10 & 
        is.character(fileLines[[1]]) & 
        {is.character(fileLines[[2]]) | is.numeric(fileLines[[2]])} & 
        is.numeric(fileLines[[3]]) & 
        is.character(fileLines[[4]]) & all(unique(fileLines[[4]]) %in% c('+', '-')) &
        {is.character(fileLines[[5]]) | is.numeric(fileLines[[5]])} & 
        is.numeric(fileLines[[6]]) & 
        is.character(fileLines[[7]]) & all(unique(fileLines[[7]]) %in% c('+', '-')) &
        is.numeric(fileLines[[8]]) & 
        is.character(fileLines[[9]]) & 
        is.character(fileLines[[10]])
    }, error = function(e) FALSE)) return('hicpro')

    # -- Check if file is Juicer format 
    # <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <cigar1> <sequence1> <mapq2> <cigar2> <sequence2> <readname1> <readname2>
    if (tryCatch({
        ncol(fileLines) >= 8 & 
        is.character(fileLines[[1]]) & sum(fileLines[[1]] == '0') > 0 & 
        {is.character(fileLines[[2]]) | is.numeric(fileLines[[2]])} & 
        is.numeric(fileLines[[3]]) & 
        is.numeric(fileLines[[4]]) & 
        is.character(fileLines[[5]]) & sum(fileLines[[5]] == '0') > 0 & 
        {is.character(fileLines[[6]]) | is.numeric(fileLines[[6]])} & 
        is.numeric(fileLines[[7]]) & 
        is.numeric(fileLines[[8]])
    }, error = function(e) FALSE)) return('juicer')

    # -- Check if file is hicstuff format 
    # <chr1> <pos1> <end1> <chr2> <start2> <end2> <frag1> <frag2> <kernel_id> <iteration> <score> <pvalue> <qvalue>
    if (tryCatch({
        ncol(fileLines) == 13 & 
        {is.character(fileLines[[1]]) | is.numeric(fileLines[[1]])} & 
        is.numeric(fileLines[[2]]) & 
        is.numeric(fileLines[[3]]) & 
        {is.character(fileLines[[4]]) | is.numeric(fileLines[[4]])} & 
        is.numeric(fileLines[[5]]) & 
        is.numeric(fileLines[[6]]) & 
        is.numeric(fileLines[[7]]) & 
        is.numeric(fileLines[[8]]) & 
        is.numeric(fileLines[[9]]) & 
        is.numeric(fileLines[[10]]) & 
        is.numeric(fileLines[[11]]) & 
        is.numeric(fileLines[[12]]) & 
        is.numeric(fileLines[[13]])
    }, error = function(e) FALSE)) return('hicstuff')

    # -- Unspecified format
    stop("Pairs file format could not be guessed. Please provide indications about the format.")
    FALSE
}