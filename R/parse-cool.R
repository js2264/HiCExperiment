#' Parsing (m)cool files
#' 
#' These functions are the workhorse internal functions used to import 
#' a `.(m)cool` file as GInteractions (wrapped into a `HiCExperiment` object
#' by `HiCExperiment()` function).
#'
#' @name parse-cool
#' 
#' @param file path to a Hi-C contact file (in (m)cool format)
#' @param resolution resolution of the contact matrix
#' @param balanced import balancing scores
#' @param pair Genomic coordinates to extract contacts for, stored as a Pairs 
#' of GRanges (e.g. 
#' S4Vectors::Pairs(GRanges("II:200000-300000"), GRanges("II:70000-100000"))). 
#' @param coords Genomic coordinates to extract contacts for, stored as a 
#' GRanges object
#' @param anchors anchors from .getCoolAnchors()
#' @param path Internal path of the cool file to check 
#' @param idx Index to extract from the cool (HDF5) file
#' @param ... Other arguments passed to .fetchCool
#' @param verbose Print resolutions in the console
#'
#' @import rhdf5
#' @import methods
#' @import InteractionSet
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr bind_rows
#' @importFrom dplyr arrange
#' @importFrom dplyr rename
#' @importFrom dplyr all_of
#' @importFrom dplyr left_join
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom IRanges subsetByOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors subjectHits
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges resize
#' @rdname parse-cool
#' @keywords internal
NULL

#' @rdname parse-cool

.getCoolAnchors <- function(file, resolution = NULL, balanced = "cooler") {
    bins <- .fetchCool(file, "bins", resolution)
    anchors <- GenomicRanges::GRanges(
        bins$chr,
        IRanges::IRanges(bins$start + 1, bins$end),
        bin_id = seq_along(bins$chr) - 1, 
        seqinfo = .cool2seqinfo(file, resolution)
    )
    names(anchors) <- paste(
        GenomicRanges::seqnames(anchors), 
        GenomicRanges::start(anchors), 
        GenomicRanges::end(anchors), sep = "_"
    )
    if ("weight" %in% names(bins) & {
        balanced == "cooler" | balanced == TRUE
    }) {
        anchors$weight <- as.numeric(bins$weight)
    } else {
        weight <- 1
    }
    return(anchors)
}

#' @rdname parse-cool

.getCountsFromPair <- function(file,
    pair,
    anchors,
    resolution = NULL
) {
    
    `%within%` <- IRanges::`%within%`

    # Make sure pair is sorted (first is first, second is after)
    pair <- sortPairs(pair)

    # Make sure pair is squared (each GRanges has same width)
    # .is_square(pair)

    # For each pair, get the coords for first and second.
    coords <- unlist(S4Vectors::zipup(pair))
    coords_list <- splitCoords(coords)

    ## Check that queried chr. exists
    if (any(!coords_list$chr %in% as.vector(GenomicRanges::seqnames(anchors)) & !is.na(coords_list$chr))) {
        sn <- paste0(unique(as.vector(GenomicRanges::seqnames(anchors))), collapse = ", ")
        stop(paste0("Some chr. are not available. Available seqnames:", sn))
    }

    ## Find out which chunks of the mcool to recover
    gr_1 <- coords[1]
    sub_1 <- which(anchors %within% gr_1)
    bin_idx_1 <- .fetchCool(file, path = "indexes/bin1_offset", resolution, idx = sub_1)
    max_bin <- max(.fetchCool(file, path = "indexes/bin1_offset", resolution))
    if ({max(bin_idx_1)+1} > max_bin) {
        chunks_1 <- seq(min(bin_idx_1)+1, max(bin_idx_1), by = 1)
    } else {
        chunks_1 <- seq(min(bin_idx_1)+1, max(bin_idx_1)+1, by = 1)
    }

    gr_2 <- coords[2]
    sub_2 <- which(anchors %within% gr_2)
    bin_idx_2 <- .fetchCool(file, path = "indexes/bin1_offset", resolution, idx = sub_2)
    max_bin <- max(.fetchCool(file, path = "indexes/bin1_offset", resolution))
    if ({max(bin_idx_2)+1} > max_bin) {
        chunks_2 <- seq(min(bin_idx_2)+1, max(bin_idx_2), by = 1)
    } else {
        chunks_2 <- seq(min(bin_idx_2)+1, max(bin_idx_2)+1, by = 1)
    }
    valid_bin2 <- unique(.fetchCool(file, path = "pixels/bin1_id", resolution, idx = chunks_2))
    
    ## Reading the chunks from the cool file
    df <- data.frame(
        bin1_id = .fetchCool(file, path = "pixels/bin1_id", resolution, idx = chunks_1),
        bin2_id = .fetchCool(file, path = "pixels/bin2_id", resolution, idx = chunks_1),
        count = .fetchCool(file, path = "pixels/count", resolution, idx = chunks_1)
    )
    
    ## Filter to only get interesting bins
    df <- df[df$bin2_id %in% valid_bin2, ]

    return(df)
}

#' @rdname parse-cool

.getCounts <- function(file,
    coords,
    anchors,
    resolution = NULL
) {
    
    `%within%` <- IRanges::`%within%`
    sub_within <- NULL
    
    ## Process coordinates
    .v <- splitCoords(coords)
    coords_chr <- .v[[1]]
    coords_start <- .v[[2]]
    coords_end <- .v[[3]]
    
    # If only chr. names are provided, find their start and stop
    if (any(is.na(coords_start) & !is.na(coords_chr))) {
        coords_start <- rep(1, length(coords_start))
        coords_end <- GenomeInfoDb::seqlengths(anchors)[coords_chr]
    }

    ## Check that queried chr. exist
    if (any(!coords_chr %in% as.vector(GenomicRanges::seqnames(anchors)) & !is.na(coords_chr))) {
        sn <- paste0(unique(as.vector(GenomicRanges::seqnames(anchors))), collapse = ", ")
        stop(paste0(
            coords_chr,
            " not in file. Available seqnames: ", 
            sn
        ))
    }

    ## Find out which chunks of the mcool to recover
    if (is.na(coords_chr)) {
        chunks <- NULL
    } else {
        gr <- GRanges(coords_chr, IRanges::IRanges(coords_start, coords_end))
        sub_within <- which(anchors %within% gr)
        bin_idx <- .fetchCool(file, path = "indexes/bin1_offset", resolution, idx = sub_within)
        max_bin <- max(.fetchCool(file, path = "indexes/bin1_offset", resolution), na.rm = TRUE)
        if ({max(bin_idx)+1} > max_bin) {
            chunks <- seq(min(bin_idx)+1, max(bin_idx), by = 1)
        } else {
            chunks <- seq(min(bin_idx)+1, max(bin_idx)+1, by = 1)
        }
    }

    ## Reading the chunks from the cool file
    df <- data.frame(
        bin1_id = .fetchCool(file, path = "pixels/bin1_id", resolution, idx = chunks),
        bin2_id = .fetchCool(file, path = "pixels/bin2_id", resolution, idx = chunks),
        count = .fetchCool(file, path = "pixels/count", resolution, idx = chunks)
    )

    ## Filter to only get interesting bins
    if (!is.null(sub_within)) {
        df <- df[df$bin2_id %in% {sub_within-1}, ]
    }

    return(df)

}

#' @rdname parse-cool

.fetchCool <- function(file, path, resolution = NULL, idx = NULL, ...) {
    path <- ifelse(
        is.null(resolution), 
        paste0("/", path), 
        paste0("/resolutions/", resolution, "/", path)
    )
    as.vector(rhdf5::h5read(
        file, name = path, index = list(idx), ..., 
        bit64conversion = 'double'
    ))
}

#' @rdname parse-cool

.dumpCool <- function(file, resolution = NULL) {
    an <- .getCoolAnchors(file, resolution)
    path <- ifelse(
        is.null(resolution), 
        '/', 
        paste0("/resolutions/", resolution)
    )
    lres <- as.vector(rhdf5::h5read(file, name = path, bit64conversion = 'double'))
    res <- list(
        'bins' = as.data.frame(lres[['bins']]),
        'chroms' = as.data.frame(lres[['chroms']]),
        'indexes' = list(
            bin1_offset = lres[['indexes']]$bin1_offset, 
            chrom_offset = lres[['indexes']]$chrom_offset
        ),
        'pixels' = as.data.frame(lres[['pixels']])
    )
    res$pixels$bin1_weight <- left_join(
       res$pixels, data.frame(an), by = c(bin1_id = 'bin_id')
    )$weight
    res$pixels$bin2_weight <- left_join(
       res$pixels, data.frame(an), by = c(bin2_id = 'bin_id')
    )$weight
    res$pixels$score <- res$pixels$count * res$pixels$bin1_weight * res$pixels$bin2_weight

    return(res)
}

#' @rdname parse-cool

.lsCoolFiles <- function(file, verbose = FALSE) {
    x <- rhdf5::h5ls(file) |> 
        dplyr::mutate(path = paste0(group, "/", name)) |> 
        dplyr::pull(path) |> 
        unique() |> 
        (\(x){gsub('//', '/', x)})()
    len <- length(x)
    if (verbose) {
        if (len > 10) {
            mess <- c(
                paste0(x[seq_len(5)], "\n"), 
                paste0("... (", len-10, " more paths)\n"), 
                paste0(x[(len-5+1):len], "\n")
            )
            message(mess)
        }
        else {
            message(x)
        }
    }
    invisible(x)
}

#' @rdname parse-cool
#' @return Silently, a numerical vector of resolutions stored in the cool file

.lsCoolResolutions <- function(file, verbose = FALSE) {
    if (.is_cool(file)) {
        x <- rhdf5::h5ls(file)
        bin_ends <- rhdf5::h5read(file, name = '/bins/end', index = list(c(1, 2)))
        res <- bin_ends[2] - bin_ends[1]
    }
    if (.is_mcool(file)) {
        x <- rhdf5::h5ls(file)
        res <- gsub("/resolutions/", "", x$group) |>
            grep(pattern = "/", invert = TRUE, value = TRUE) |>
            unique() |>
            as.numeric() |>
            sort() |>
            as.character()
    }
    if (verbose) message(S4Vectors::coolcat("resolutions(%d): %s", res))
    invisible(as.integer(res))
}

#' @rdname parse-cool

.cool2seqinfo <- function(file, resolution = NULL) {
    chroms <- .fetchCool(file, "chroms", resolution)
    seqinfo <- GenomeInfoDb::Seqinfo(
        seqnames = as.vector(chroms$name),
        seqlengths = as.vector(chroms$length)
    )
    return(seqinfo)
}

#' @rdname parse-cool

.cool2gi <- function(file, coords = NULL, resolution = NULL) {

    # Mutate Pairs provided as characters to real Pairs
    if (!is.null(coords)) {
        if (grepl('\\|', coords)) {
            si <- .cool2seqinfo(file, resolution)
            if (
                grepl(
                    '[A-Za-z0-9]*\\|[A-Za-z0-9]*$', coords
                )
            ) { # e.g. 'II|III'
                chr1 <- strsplit(coords, '\\|')[[1]][1]
                chr2 <- strsplit(coords, '\\|')[[1]][2]
                if (!all(c(chr1, chr2) %in% seqnames(si))) {
                    stop("One or all of the provided seqnames is not found.")
                }
                gr1 <- as(si[chr1], 'GRanges')
                GenomeInfoDb::seqlevels(gr1) <- GenomeInfoDb::seqlevels(si)
                gr2 <- as(si[chr2], 'GRanges')
                GenomeInfoDb::seqlevels(gr1) <- GenomeInfoDb::seqlevels(si)
                coords <- S4Vectors::Pairs(
                    sort(c(gr1, gr2))[1], 
                    sort(c(gr1, gr2))[2]
                )
            }
            else {
                coords <- char2coords(coords)
                gr1 <- S4Vectors::first(coords)
                GenomeInfoDb::seqlevels(gr1) <- GenomeInfoDb::seqlevels(si)
                gr2 <- S4Vectors::second(coords)
                GenomeInfoDb::seqlevels(gr1) <- GenomeInfoDb::seqlevels(si)
                coords <- S4Vectors::Pairs(
                    sort(c(gr1, gr2))[1], 
                    sort(c(gr1, gr2))[2]
                )
            }

        }
    }

    # Check if the provided coords are GRanges or Pairs
    is_pair <- is(coords, 'Pairs')

    # Get anchors from mcool
    anchors <- .getCoolAnchors(file, resolution)

    # Get raw counts for bins from mcool
    if (!is_pair) {
        cnts <- .getCounts(file, coords = coords, anchors = anchors, resolution = resolution)
    }
    else {
        cnts <- .getCountsFromPair(file, pair = coords, anchors = anchors, resolution = resolution)
    }

    # Associate raw counts for bins to corresponding anchors
    gi <- InteractionSet::GInteractions(
        anchors[cnts$bin1_id+1],
        anchors[cnts$bin2_id+1],
        count = cnts$count
    )
    gi$bin_id1 <- cnts$bin1_id
    gi$bin_id2 <- cnts$bin2_id
    
    # Get balanced counts if they exist
    if (!is.null(gi$anchor1.weight) & !is.null(gi$anchor2.weight)) {
        gi$score <- gi$count * gi$anchor1.weight * gi$anchor2.weight
    } 
    else {
        gi$score <- gi$count
    }

    # Fix regions by adding the empty ones with no original 'count'
    gi <- .fixRegions(gi, anchors, coords)

    return(gi)
}
