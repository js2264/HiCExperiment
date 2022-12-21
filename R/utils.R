#' @title Utils functions 
#' 
#' @name HiCExperiment utils
#' 
#' @rdname utils
#' 
#' @description 
#' 
#' Utilities to facilitate parsing/handling of coordinates, GInteractions, 
#' Pairs, ...
#'
#' @return Reformatted coordinates or GInteractions
#' 
#' @param coords coords
#' @import stringr
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end

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

#' @param big.mark big.mark
#' @rdname utils

coords2char <- function(coords, big.mark = ',') {
    if (is(coords, 'GRanges')) {
        chr <- as.vector(GenomicRanges::seqnames(coords))
        start <- GenomicRanges::start(coords)
        end <- GenomicRanges::end(coords)
        paste0(chr, ':', format(start, big.mark = big.mark), '-', format(end, big.mark = big.mark))
    }
    else {
        if (grepl('\\|', coords)) {
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

#' @param char char (e.g. "II:30000-50000" or "II:30000-50000|II:60000-80000")
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
        '[A-Za-z0-9]*:[0-9]*-[0-9]*\\|[A-Za-z0-9]*:[0-9]*-[0-9]*$', 
        char
    )) {
        splitst <- stringr::str_split(char, '\\|')[[1]]
        S4Vectors::Pairs(
            GenomicRanges::GRanges(splitst[[1]]), 
            GenomicRanges::GRanges(splitst[[2]])
        )
    }
    else if (grepl('[A-Za-z0-9]*:[0-9]*[-:][0-9]*$', char)) {
        S4Vectors::Pairs(
            GenomicRanges::GRanges(char), 
            GenomicRanges::GRanges(char)
        )
    }
    else {
        stop("Cannot coerce string into a Pairs object")
    }
}

#' @param pairs pairs
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

#' @param df df
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

#' @param A matrix
#' @param k secondary diagonal k
#' @rdname utils

`sdiag` <- function(A, k = 0) {
    p <- ncol(A)
    n <- nrow(A)
    if (k>p-1||-k > n-1) return()
    if (k >= 0) {
        i <- seq_len(n)
        j <- (k+1):p
    } 
    else {
        i <- (-k+1):n
        j <- seq_len(p)
    }
    if (length(i)>length(j)) i <- i[seq_along(j)] else j <- j[seq_along(i)]
    ii <- i + (j-1) * n 
    A[ii]
}

#' @param diag vector of distances to diagonal
#' @param score scores to parse into symmetrical matrix
#' @rdname utils

.df_to_symmmat <- function(diag, score) {
    .n <- length(diag)
    .l <- vector(mode = 'list', length = .n)
    .l[[1]] <- score
    .l[-1] <- lapply(seq(2, .n), function(k) {c(tail(score, k-1), head(score, 101-k+1))})
    .l <- do.call(rbind, .l)
    .l[lower.tri(.l)] <- .l[upper.tri(.l)]
    return(.l)
}

#' @param dump dumped contacts, e.g. from .dumpCool
#' @param threshold maximum distance to compute distance decay for
#' @rdname utils

distance_decay <- function(dump, threshold = 1e10) {
    resolution <- as.vector(dump$bins[1, ]$end - dump$bins[1, ]$start)
    dump$bins$bin_id <- seq(0, nrow(dump$bins)-1)
    df <- dump$pixels
    j1 <- dplyr::left_join(
        dump$pixels, dump$bins, by = c(bin1_id = 'bin_id')
    )
    j2 <- dplyr::left_join(
        dump$pixels, dump$bins, by = c(bin2_id = 'bin_id')
    )
    df$chrom1 <- j1$chrom
    df$weight1 <- j1$weight
    df$chrom2 <- j2$chrom
    df$weight2 <- j2$weight
    df <- dplyr::filter(df, chrom1 == chrom2)
    df <- dplyr::mutate(df, 
        diag = {as.vector(bin2_id) - as.vector(bin1_id)}, 
        distance = diag * resolution
    ) |> 
        dplyr::filter(distance <= threshold)
    mod <- dplyr::group_by(df, diag, distance) |> 
        dplyr::mutate(score = count * weight1 * weight2) |> 
        dplyr::summarize(score = mean(score, na.rm = TRUE))
    return(mod)
}
