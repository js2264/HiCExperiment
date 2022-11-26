#' Utils functions 
#' 
#' Utilities to facilitate parsing/handling of coordinates, GInteractions, 
#' Pairs, ...
#'
#' @param coords coords
#' @return a list containing `chr`, `start` and `end`
#' 
#' @import stringr
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @rdname utils

splitCoords <- function(coords) {
    if (is(coords, 'GRanges')) {
        chr <- as.vector(GenomicRanges::seqnames(coords))
        start <- GenomicRanges::start(coords)
        end <- GenomicRanges::end(coords)
        list(
            "chr" = chr,
            "start" = start,
            "end" = end
        )
    }
    else {
        chr <- stringr::str_replace(coords, ":.*", "")
        chr <- ifelse(length(chr) == 0, NA, chr)
        start <- suppressWarnings(as.numeric(stringr::str_replace_all(coords, ".*:|-.*", "")))
        start <- ifelse(length(start) == 0, NA, start)
        end <- suppressWarnings(as.numeric(stringr::str_replace(coords, ".*-", "")))
        end <- ifelse(length(end) == 0, NA, end)
        list(
            "chr" = chr,
            "start" = start,
            "end" = end
        )
    }
}

#' coords2char
#'
#' @param coords coords
#' @param big.mark big.mark
#' @return a character string
#'
#' @import stringr
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @rdname utils

coords2char <- function(coords, big.mark = ',') {
    if (is(coords, 'GRanges')) {
        chr <- as.vector(GenomicRanges::seqnames(coords))
        start <- GenomicRanges::start(coords)
        end <- GenomicRanges::end(coords)
        paste0(chr, ':', format(start, big.mark = big.mark), '-', format(end, big.mark = big.mark))
    }
    else {
        if (grepl('x', coords)) {
            coords
        }
        else {
            chr <- stringr::str_replace(coords, ":.*", "")
            start <- suppressWarnings(as.numeric(stringr::str_replace_all(coords, ".*:|-.*", "")))
            end <- suppressWarnings(as.numeric(stringr::str_replace(coords, ".*-", "")))
            if (is.na(start)) {
                return(chr)
            }
            paste0(chr, ':', format(start, big.mark = big.mark, scientific = FALSE), '-', format(end, big.mark = big.mark, scientific = FALSE))
        }
    }
}

#' char2coords
#'
#' @param char char (e.g. "II:30000-50000" or "II:30000-50000 x II:60000-80000")
#' @return a S4Vectors::Pairs object
#'
#' @import stringr
#' @importFrom methods is
#' @importFrom stringr str_split
#' @importFrom S4Vectors Pairs
#' @importFrom GenomicRanges GRanges
#' @rdname utils

char2coords <- function(char) {
    # if (methods::is(char, 'Pairs')) {
    #     return(char)
    # }
    if (grepl(
        '[A-Za-z0-9]*:[0-9]*-[0-9]* [xX/-;] [A-Za-z0-9]*:[0-9]*-[0-9]*$', 
        char
    )) {
        splitst <- stringr::str_split(char, ' . ')[[1]]
        S4Vectors::Pairs(
            GenomicRanges::GRanges(splitst[[1]]), 
            GenomicRanges::GRanges(splitst[[2]])
        )
    }
    else if (grepl('[A-Za-z0-9]*:[0-9]*-[0-9]*$', char)) {
        S4Vectors::Pairs(
            GenomicRanges::GRanges(char), 
            GenomicRanges::GRanges(char)
        )
    }
    else {
        stop("Cannot coerce string into a Pairs object")
    }
}

#' fullContactInteractions
#'
#' @param chr chr
#' @param start start
#' @param end end
#' @param binning binning
#' @return a GenomicInteractions object
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom InteractionSet GInteractions
#' @rdname utils

fullContactInteractions <- function(chr, start, end, binning) {
    full_anchors <- GenomicRanges::GRanges(
        seqnames = chr, 
        IRanges::IRanges(
            start = seq(start, end-1, by = binning),
            width = binning
        )
    )
    InteractionSet::GInteractions(
        full_anchors[rep(seq_along(full_anchors), length(full_anchors))], 
        full_anchors[rep(seq_along(full_anchors), each = length(full_anchors))], 
        full_anchors
    )
}

#' sortPairs
#'
#' @param pairs pairs
#' @return a Pairs object
#'
#' @importFrom S4Vectors zipup
#' @importFrom S4Vectors zipdown
#' @importFrom GenomicRanges GRangesList
#' @rdname utils

sortPairs <- function(pairs) {
    p <- S4Vectors::zipup(pairs)
    p_sorted <- GenomicRanges::GRangesList(lapply(p, function(gr) {
        sort(gr)
    }))
    S4Vectors::zipdown(p_sorted)
}

#' asGInteractions
#'
#' @param df df
#' @return a GenomicInteractions object
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom InteractionSet GInteractions
#' @rdname utils

asGInteractions <- function(df) {
    gi <- InteractionSet::GInteractions(
        anchor1 = GenomicRanges::GRanges(
            df$seqnames1, IRanges::IRanges(df$start1, df$end1)
        ),
        anchor2 = GenomicRanges::GRanges(
            df$seqnames2, IRanges::IRanges(df$start2, df$end2)
        )
    )
    if ('score' %in% colnames(df)) gi$score <- df$score
    gi
}
