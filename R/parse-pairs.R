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

pairs2gi <- function(
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
        file_lines <- vroom::vroom(
            file,
            n_max = 100000,
            col_names = FALSE,
            comment = '#',
            progress = FALSE, 
            show_col_types = FALSE
        )
        chr1.field <- 2
        start1.field <- 3
        chr2.field <- 4
        start2.field <- 5
        strand1.field <- 6
        strand2.field <- 7
        if (ncol(file_lines) >= 9) {
            frag1.field <- 8
            frag2.field <- 9
        }
    }
    if (.fmt == 'hicpro') {
        # <ID> <chr1> <pos1> <str1> <chr2> <pos2> <str2> <isize> <frag1> <frag2> <mapq1> <mapq2> [<allele-spe. info>]
        file_lines <- vroom::vroom(
            file,
            n_max = 100000,
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
        if (ncol(file_lines) >= 10) {
            frag1.field <- 8
            frag2.field <- 9
        }
    }
    if (any(c(
        is.null(chr1.field), 
        is.null(start1.field), 
        is.null(strand1.field), 
        is.null(chr2.field), 
        is.null(start2.field), 
        is.null(strand2.field)
    ))) stop("Impossible to parse pairs file")

    # -- Import pairs
    sel1 <- (c(chr1.field, start1.field, strand1.field, frag1.field))
    sel2 <- (c(chr2.field, start2.field, strand2.field, frag2.field))
    anchors1 <- vroom::vroom(
        file,
        n_max = nrows,
        col_select = dplyr::all_of(sel1),
        comment = '#',
        col_names = FALSE,
        show_col_types = FALSE, 
        num_threads = nThread
    )
    anchors2 <- vroom::vroom(
        file,
        n_max = nrows,
        col_select = dplyr::all_of(sel2),
        comment = '#',
        col_names = FALSE,
        show_col_types = FALSE, 
        num_threads = nThread
    )  
    anchor_one <- GenomicRanges::GRanges(
        anchors1[[1]],
        IRanges::IRanges(anchors1[[2]], width = 1), 
        strand = anchors1[[3]]
    )
    anchor_two <- GenomicRanges::GRanges(
        anchors2[[1]],
        IRanges::IRanges(anchors2[[2]], width = 1), 
        strand = anchors2[[3]]
    )
    gi <- InteractionSet::GInteractions(anchor_one, anchor_two)
    
    if (!is.null(frag1.field) & !is.null(frag2.field)) {
        gi$frag1 <- anchors1[[4]]
        gi$frag2 <- anchors2[[4]]
    } 
    else {
        gi$frag1 <- NA
        gi$frag2 <- NA
    }
    gi$distance <- InteractionSet::pairdist(gi) 
    return(gi)
}

.guess_pairs_format <- function(file) {

    lines <- readLines(file, n = 1000)
    file_comments <- grep('^#', lines, value = TRUE)
    file_lines <- vroom::vroom(
        file,
        n_max = 100000,
        col_names = FALSE,
        comment = '#',
        progress = FALSE, 
        show_col_types = FALSE
    )
    ncols <- ncol(file_lines)

    # -- Check if file is .pairs format (4DN)
    # <ID> <chr1> <pos1> <chr2> <pos2> <str1> <str2> [<frag1> <frag2>]
    if (any(grepl('pairs format', file_comments))) return('4dn')
    if (tryCatch({
        ncol(file_lines) >= 7 & 
        is.character(file_lines[[1]]) & # <ID>
        {is.character(file_lines[[2]]) | is.numeric(file_lines[[2]])} & # <chr1> 
        is.numeric(file_lines[[3]]) & # <pos1> 
        {is.character(file_lines[[4]]) | is.numeric(file_lines[[4]])} & # <chr2> 
        is.numeric(file_lines[[5]]) & # <pos2> 
        is.character(file_lines[[6]]) & all(unique(file_lines[[6]]) %in% c('+', '-')) & # <str1> 
        is.character(file_lines[[7]]) & all(unique(file_lines[[7]]) %in% c('+', '-')) # <str2>
    }, error = function(e) FALSE)) return('4dn')

    # -- Check if file is .validPairs format (HiC-Pro)
    # <ID> <chr1> <pos1> <str1> <chr2> <pos2> <str2> <isize> <frag1> <frag2> <mapq1> <mapq2> [<allele-spe. info>]
    if (tryCatch({
        ncol(file_lines) >= 10 & 
        is.character(file_lines[[1]]) & 
        {is.character(file_lines[[2]]) | is.numeric(file_lines[[2]])} & 
        is.numeric(file_lines[[3]]) & 
        is.character(file_lines[[4]]) & all(unique(file_lines[[4]]) %in% c('+', '-')) &
        {is.character(file_lines[[5]]) | is.numeric(file_lines[[5]])} & 
        is.numeric(file_lines[[6]]) & 
        is.character(file_lines[[7]]) & all(unique(file_lines[[7]]) %in% c('+', '-')) &
        is.numeric(file_lines[[8]]) & 
        is.character(file_lines[[9]]) & 
        is.character(file_lines[[10]])
    }, error = function(e) FALSE)) return('hicpro')

    # -- Check if file is Juicer format 
    # <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <cigar1> <sequence1> <mapq2> <cigar2> <sequence2> <readname1> <readname2>
    if (tryCatch({
        ncol(file_lines) >= 8 & 
        is.character(file_lines[[1]]) & sum(file_lines[[1]] == '0') > 0 & 
        {is.character(file_lines[[2]]) | is.numeric(file_lines[[2]])} & 
        is.numeric(file_lines[[3]]) & 
        is.numeric(file_lines[[4]]) & 
        is.character(file_lines[[5]]) & sum(file_lines[[5]] == '0') > 0 & 
        {is.character(file_lines[[6]]) | is.numeric(file_lines[[6]])} & 
        is.numeric(file_lines[[7]]) & 
        is.numeric(file_lines[[8]])
    }, error = function(e) FALSE)) return('juicer')

    # -- Unspecified format
    stop("Pairs file format could not be guessed. Please provide indications about the format.")
    FALSE
}