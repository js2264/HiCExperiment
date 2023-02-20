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
#' These functions are not exported.
#'
#' @return Reformatted coordinates or GInteractions.
#' 
#' @param coords A set of genomic coordinates (either as a GRanges
#' object or as a character string)
#' @param big.mark Separator for thousands when printing out genomic coordinates 
#' as character
#' @param char char (e.g. "II:30000-50000" or "II:30000-50000|II:60000-80000")
#' @param pairs Pairs object
#' @param df a data.frame to turn into a GInteraction object. 
#' @param A Numerical matrix
#' @param k secondary diagonal k
#' @param diag vector of distances to diagonal
#' @param score scores to parse into symmetrical matrix
#' @param dump dumped contacts as GInteractions, e.g. from .dumpCool
#' @param threshold maximum distance to compute distance decay for
#' @param file path to a HiC contact matrix file
#' @param resolution Resolution to use with the HiC contact matrix file
#' @param gis GInteractions object
#' @param bins Larger set of regions (usually bins from HiCExperiment)
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @keywords internal
NULL

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
        # chr <- stringr::str_replace(coords, ":.*", "")
        chr <- gsub(":.*", "", coords)
        chr <- ifelse(length(chr) == 0, NA, chr)
        # start <- suppressWarnings(as.numeric(stringr::str_replace_all(coords, ".*:|-.*", "")))
        start <- suppressWarnings(as.numeric((gsub(".*:|-.*", "", coords))))
        start <- ifelse(length(start) == 0, NA, start)
        # end <- suppressWarnings(as.numeric(stringr::str_replace(coords, ".*-", "")))
        end <- suppressWarnings(as.numeric(gsub(".*-", "", coords)))
        end <- ifelse(length(end) == 0, NA, end)
        list(
            "chr" = chr,
            "start" = start,
            "end" = end
        )
    }
}

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
            # chr <- stringr::str_replace(coords, ":.*", "")
            chr <- gsub(":.*", "", coords)
            # start <- suppressWarnings(as.numeric(stringr::str_replace_all(coords, ".*:|-.*", "")))
            start <- suppressWarnings(as.numeric((gsub(".*:|-.*", "", coords))))
            # end <- suppressWarnings(as.numeric(stringr::str_replace(coords, ".*-", "")))
            end <- suppressWarnings(as.numeric(gsub(".*-", "", coords)))
            if (is.na(start)) {
                return(chr)
            }
            paste0(chr, ':', format(start, big.mark = big.mark, scientific = FALSE), '-', format(end, big.mark = big.mark, scientific = FALSE))
        }
    }
}

#' @importFrom methods is
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
        splitst <- strsplit(char, '\\|')[[1]]
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
        stop("Cannot coerce the provided string into a Pairs object")
    }
}

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
    if ('bin_id1' %in% colnames(df) & 'bin_id2' %in% colnames(df)) {
        gi$bin_id1 <- df$bin_id1
        gi$bin_id2 <- df$bin_id2
    }
    gi
}

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

#' @rdname utils

.df2symmmat <- function(diag, score) {
    .n <- length(diag)
    .l <- vector(mode = 'list', length = .n)
    .l[[1]] <- score
    .l[-1] <- lapply(seq(2, .n), function(k) {
        c(utils::tail(score, k-1), utils::head(score, .n-k+1))
    })
    .l <- do.call(rbind, .l)
    .l[lower.tri(.l)] <- t(.l)[lower.tri(.l)]
    return(.l)
}

#' @rdname utils

distanceDecay <- function(dump, threshold = NULL) {
    resolution <- as.vector(dump$bins[1, ]$end - dump$bins[1, ]$start)
    dump$bins$bin_id <- seq(0, nrow(dump$bins)-1)
    df <- dump$pixels
    if (!'score' %in% colnames(df)) {
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
    }
    df <- dplyr::mutate(df, 
        diag = {as.vector(bin2_id) - as.vector(bin1_id)}, 
        distance = diag * resolution
    )
    if (!is.null(threshold)) df <- dplyr::filter(df, distance <= threshold)
    mod <- dplyr::group_by(df, diag, distance) |> 
        dplyr::summarize(score = mean(score, na.rm = TRUE))
    return(mod)
}

#' @rdname utils

detrendingModel <- function(file, resolution) {
    if (.is_cool(file)) {
        l <- .dumpCool(file, resolution = NULL)
    }
    if (.is_mcool(file)) {
        l <- .dumpCool(file, resolution = resolution)
    }
    else if (.is_hic(file)) {
        l <- .dumpHic(file, resolution = resolution)
    }
    else if (.is_hicpro_matrix(file)) {
        l <- .dumpHicpro(file, bed)
    }
    detrendingModel <- distanceDecay(l)
    return(detrendingModel)
}

#' @rdname utils

.fixRegions <- function(gis, bins, coords) {

    # Different possible situations: 
    #   NULL
    #   'II:10000-20000'
    #   'II'
    #   'II:10000-20000|III:50000-90000'
    #   'II|III'
    #   c('II', 'III')

    if (is.null(coords)) { ## genome-wide 
        re <- bins
    }
    else if (
        is(coords, 'Pairs')
    ) { ## coords already in Pairs 
        valid_bins <- c(
            subsetByOverlaps(bins, S4Vectors::first(coords))$bin_id, 
            subsetByOverlaps(bins, S4Vectors::second(coords))$bin_id
        )
        re <- bins[bins$bin_id %in% valid_bins]
    }
    else if (
        length(coords) > 1
    ) { # c('II', 'III')
        valid_bins <- as.vector(seqnames(bins)) %in% coords
        re <- bins[valid_bins]
    }
    else if (
        grepl(
            '[A-Za-z0-9]*:[0-9]*-[0-9]*\\|[A-Za-z0-9]*:[0-9]*-[0-9]*$', coords
        ) | grepl(
            '[A-Za-z0-9]*:[0-9]*-[0-9]*$', coords
        ) 
    ) { # e.g. 'II:10000-20000' or 'II:10000-20000|III:50000-90000'
        coords <- char2coords(coords)
        valid_bins <- c(
            subsetByOverlaps(bins, S4Vectors::first(coords))$bin_id, 
            subsetByOverlaps(bins, S4Vectors::second(coords))$bin_id
        )
        re <- bins[bins$bin_id %in% valid_bins]
    }
    else if (
        grepl(
            '[A-Za-z0-9]*\\|[A-Za-z0-9]*$', coords
        )
    ) { # e.g. 'II|III'
        chr1 <- strsplit(coords, '\\|')[[1]][1]
        chr2 <- strsplit(coords, '\\|')[[1]][2]
        if (!all(c(chr1, chr2) %in% seqnames(GenomeInfoDb::seqinfo(gis)))) {
            stop("One or all of the provided seqnames is not found.")
        }
        si <- GenomeInfoDb::seqinfo(gis)
        gr1 <- as(si[chr1], 'GRanges')
        gr2 <- as(si[chr2], 'GRanges')
        valid_bins <- c(
            subsetByOverlaps(bins, gr1)$bin_id, 
            subsetByOverlaps(bins, gr2)$bin_id
        )
        re <- bins[bins$bin_id %in% valid_bins]
    }
    else if (
        all(coords %in% seqnames(GenomeInfoDb::seqinfo(gis)))
    ) { # e.g. 'II'
        valid_bins <- bins$bin_id[as.vector(seqnames(bins)) %in% coords]
        re <- bins[bins$bin_id %in% valid_bins]
    }

    re$chr <- GenomicRanges::seqnames(re)
    re$center <- GenomicRanges::start(GenomicRanges::resize(re, fix = "center", width = 1))
    re$bin_id <- bins$bin_id[BiocGenerics::match(re, bins)]
    an <- anchors(gis)
    gis <- gis[IRanges::overlapsAny(an[[1]], re) & IRanges::overlapsAny(an[[2]], re)]
    InteractionSet::replaceRegions(gis) <- re
    return(gis)

}
