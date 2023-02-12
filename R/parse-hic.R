#' Parsing hic files
#' 
#' These functions are the workhorse internal functions used to import 
#' a `.hic` file as GInteractions (wrapped into a `HiCExperiment` object
#' by `HiCExperiment()` function).
#'
#' @param file path to a Hi-C contact file in .hic format
#' @param coords NULL, character, or GRanges. 
#'   Can also be a Pairs object of paired GRanges (length of 1).
#' @param resolution resolution of the contact matrix to use
#' @return a GInteractions object
#'
#' @import InteractionSet
#' @import strawr
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges resize
#' @rdname parse-hic

.hic2gi <- function(file, coords = NULL, resolution = NULL) {
    
    file <- gsub('~', Sys.getenv('HOME'), file)

    # Mutate Pairs provided as characters to real Pairs
    if (!is.null(coords)) {
        if (grepl('\\|', coords)) {
            si <- .hic2seqinfo(file)
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

    # Get anchors from hic
    anchors <- .getHicAnchors(file, resolution)
    si <- GenomeInfoDb::seqinfo(anchors)

    # Get raw counts for bins from hic
    if (is.null(coords)) {
        combs <- combn(GenomeInfoDb::seqlevels(anchors), m = 2) |> 
                as.data.frame() |> 
                t()
        colnames(combs) <- c('one', 'two')
        combs <- rbind(
            data.frame(
                'one' = GenomeInfoDb::seqlevels(anchors), 
                'two' = GenomeInfoDb::seqlevels(anchors)
            ),
            combs
        )
        gis_ <- apply(combs, 1, function(coords) {tryCatch(expr = {
            parsed_hic <- strawr::straw(
                fname = file, 
                binsize = resolution, 
                chr1loc = coords[[1]], 
                chr2loc = coords[[2]], 
                unit = 'BP', 
                matrix = "observed", 
                norm = 'NONE'
            )
            parsed_hic_balanced <- strawr::straw(
                fname = file, 
                binsize = resolution, 
                chr1loc = coords[[1]], 
                chr2loc = coords[[2]], 
                unit = 'BP', 
                matrix = "observed", 
                norm = 'KR'
            )
            seqnames1 <- gsub(':.*', '', coords[[1]])
            seqnames2 <- gsub(':.*', '', coords[[2]])
            an1 <- GenomicRanges::GRanges(
                seqnames = seqnames1, 
                ranges = IRanges::IRanges(
                    start = parsed_hic$x + 1, 
                    width = resolution
                )
            )
            ean1 <- GenomicRanges::end(an1)
            maxean1 <- GenomeInfoDb::seqlengths(si)[GenomeInfoDb::seqlevels(si) == coords[[1]]]
            GenomicRanges::end(an1)[ean1 > maxean1] <- maxean1
            GenomeInfoDb::seqlevels(an1) <- GenomeInfoDb::seqlevels(si)
            GenomeInfoDb::seqinfo(an1) <- si
            an2 <- GenomicRanges::GRanges(
                seqnames = seqnames2, 
                ranges = IRanges::IRanges(
                    start = parsed_hic$y + 1, 
                    width = resolution
                )
            )
            ean2 <- GenomicRanges::end(an2)
            maxean2 <- GenomeInfoDb::seqlengths(si)[GenomeInfoDb::seqlevels(si) == coords[[2]]]
            GenomicRanges::end(an2)[ean2 > maxean2] <- maxean2
            GenomeInfoDb::seqlevels(an2) <- GenomeInfoDb::seqlevels(si)
            GenomeInfoDb::seqinfo(an2) <- si
            gi <- InteractionSet::GInteractions(
                an1, 
                an2, 
            )

            # Associate raw counts for bins to corresponding anchors
            gi$count <- parsed_hic$counts
            gi$score <- parsed_hic_balanced$counts

            # Return gi
            as.data.frame(gi)
        }, error = function(e) {
            as.data.frame(InteractionSet::GInteractions())
        })})
        gis_ <- gis_[sapply(gis_, nrow) > 0]
        full_parsed_hic <- do.call(rbind, gis_)
        an1 <- GenomicRanges::GRanges(
            seqnames = full_parsed_hic$seqnames1, 
            ranges = IRanges::IRanges(
                start = full_parsed_hic$start1, 
                width = resolution
            )
        )
        suppressWarnings(GenomeInfoDb::seqinfo(an1) <- si)
        an1 <- GenomicRanges::trim(an1)
        an2 <- GenomicRanges::GRanges(
            seqnames = full_parsed_hic$seqnames2, 
            ranges = IRanges::IRanges(
                start = full_parsed_hic$start2, 
                width = resolution
            )
        )
        suppressWarnings(GenomeInfoDb::seqinfo(an2) <- si)
        an2 <- GenomicRanges::trim(an2)
        gi <- InteractionSet::GInteractions(
            an1, 
            an2, 
        )

        # Associate raw counts for bins to corresponding anchors
        gi$count <- full_parsed_hic$count
        gi$score <- full_parsed_hic$score
    }
    else if (!is_pair) {
        parsed_hic <- strawr::straw(
            fname = file, 
            binsize = resolution, 
            chr1loc = coords, 
            chr2loc = coords, 
            unit = 'BP', 
            matrix = "observed", 
            norm = 'NONE'
        )
        parsed_hic_balanced <- strawr::straw(
            fname = file, 
            binsize = resolution, 
            chr1loc = coords, 
            chr2loc = coords, 
            unit = 'BP', 
            matrix = "observed", 
            norm = 'KR'
        )
        seqnames <- gsub(':.*', '', coords)
        an1 <- GenomicRanges::GRanges(
            seqnames = seqnames, 
            ranges = IRanges::IRanges(
                start = parsed_hic$x + 1, 
                width = resolution
            )
        )
        ean1 <- GenomicRanges::end(an1)
        maxean1 <- GenomeInfoDb::seqlengths(si)[GenomeInfoDb::seqlevels(si) == seqnames]
        GenomicRanges::end(an1)[ean1 > maxean1] <- maxean1
        GenomeInfoDb::seqlevels(an1) <- GenomeInfoDb::seqlevels(si)
        GenomeInfoDb::seqinfo(an1) <- si
        an2 <- GenomicRanges::GRanges(
            seqnames = seqnames, 
            ranges = IRanges::IRanges(
                start = parsed_hic$y + 1, 
                width = resolution
            )
        )
        ean2 <- GenomicRanges::end(an2)
        maxean2 <- GenomeInfoDb::seqlengths(si)[GenomeInfoDb::seqlevels(si) == seqnames]
        GenomicRanges::end(an2)[ean2 > maxean2] <- maxean2
        GenomeInfoDb::seqlevels(an2) <- GenomeInfoDb::seqlevels(si)
        GenomeInfoDb::seqinfo(an2) <- si
        gi <- InteractionSet::GInteractions(
            an1, 
            an2, 
        )

        # Associate raw counts for bins to corresponding anchors
        gi$count <- parsed_hic$counts
        gi$score <- parsed_hic_balanced$counts
    }
    else {
        parsed_hic <- strawr::straw(
            fname = file, 
            binsize = resolution, 
            chr1loc = as.character(S4Vectors::first(coords)), 
            chr2loc = as.character(S4Vectors::second(coords)), 
            unit = 'BP', 
            matrix = "observed", 
            norm = 'NONE'
        )
        parsed_hic_balanced <- strawr::straw(
            fname = file, 
            binsize = resolution, 
            chr1loc = as.character(S4Vectors::first(coords)), 
            chr2loc = as.character(S4Vectors::second(coords)), 
            unit = 'BP', 
            matrix = "observed", 
            norm = 'KR'
        )
        seqnames1 <- gsub(':.*', '', as.character(S4Vectors::first(coords)))
        seqnames2 <- gsub(':.*', '', as.character(S4Vectors::second(coords)))
        an1 <- GenomicRanges::GRanges(
            seqnames = seqnames1, 
            ranges = IRanges::IRanges(
                start = parsed_hic$x + 1, 
                width = resolution
            )
        )
        ean1 <- GenomicRanges::end(an1)
        maxean1 <- GenomeInfoDb::seqlengths(si)[GenomeInfoDb::seqlevels(si) == seqnames1]
        GenomicRanges::end(an1)[ean1 > maxean1] <- maxean1
        GenomeInfoDb::seqlevels(an1) <- GenomeInfoDb::seqlevels(si)
        GenomeInfoDb::seqinfo(an1) <- si
        an2 <- GenomicRanges::GRanges(
            seqnames = seqnames2, 
            ranges = IRanges::IRanges(
                start = parsed_hic$y + 1, 
                width = resolution
            )
        )
        ean2 <- GenomicRanges::end(an2)
        maxean2 <- GenomeInfoDb::seqlengths(si)[GenomeInfoDb::seqlevels(si) == seqnames2]
        GenomicRanges::end(an2)[ean2 > maxean2] <- maxean2
        GenomeInfoDb::seqlevels(an2) <- GenomeInfoDb::seqlevels(si)
        GenomeInfoDb::seqinfo(an2) <- si
        gi <- InteractionSet::GInteractions(
            an1, 
            an2, 
        )

        # Associate raw counts for bins to corresponding anchors
        gi$count <- parsed_hic$counts
        gi$score <- parsed_hic_balanced$counts
    }

    # Match anchor ID with each interaction
    gi$bin_id1 <- S4Vectors::subjectHits(
        GenomicRanges::findOverlaps(an1, anchors)
    ) - 1
    gi$bin_id2 <- S4Vectors::subjectHits(
        GenomicRanges::findOverlaps(an2, anchors)
    ) - 1
    
    # Fix regions by adding the empty ones with no original 'count'
    gi <- .fixRegions(gi, anchors, coords)
    
    return(gi)
}

#' @param verbose Print resolutions in the console
#' @return vector
#' @rdname parse-hic

lsHicResolutions <- function(file, verbose = FALSE) {
    res <- rev(strawr::readHicBpResolutions(file))
    if (verbose) message(S4Vectors::coolcat("resolutions(%d): %s", res))
    invisible(as.integer(res))
}

#' @rdname parse-hic

.getHicAnchors <- function(file, resolution = NULL) {
    si <- .hic2seqinfo(file)
    anchors <- GenomicRanges::tileGenome(
        si, tilewidth = resolution, cut.last.tile.in.chrom = TRUE
    )
    anchors$bin_id <- seq_along(anchors) - 1
    names(anchors) <- paste(
        GenomicRanges::seqnames(anchors), 
        GenomicRanges::start(anchors), 
        GenomicRanges::end(anchors), 
        sep = "_"
    )
    weight <- 1
    return(anchors)
}

#' @rdname parse-hic

.hic2seqinfo <- function(file) {
    strawr::readHicChroms(file) |> 
        setNames(c('seqnames', 'seqlengths')) |>
        dplyr::filter(!seqnames %in% c('ALL', 'All')) |>
        as("Seqinfo")
}

#' @rdname parse-hic

.dumpHic <- function(file, resolution = NULL) {
    check_hic_format(file, resolution)
    anchors <- .getHicAnchors(file, resolution)
    bins <- as.data.frame(anchors)
    si <- GenomeInfoDb::seqinfo(anchors)
    combs <- data.frame(
        'one' = GenomeInfoDb::seqlevels(anchors), 
        'two' = GenomeInfoDb::seqlevels(anchors)
    )
    pixs <- apply(combs, 1, function(coords) {tryCatch(expr = {
        parsed_hic <- strawr::straw(
            fname = file, 
            binsize = resolution, 
            chr1loc = coords[[1]], 
            chr2loc = coords[[2]], 
            unit = 'BP', 
            matrix = "observed", 
            norm = 'NONE'
        )
        parsed_hic_balanced <- strawr::straw(
            fname = file, 
            binsize = resolution, 
            chr1loc = coords[[1]], 
            chr2loc = coords[[2]], 
            unit = 'BP', 
            matrix = "observed", 
            norm = 'KR'
        )
        seqnames1 <- gsub(':.*', '', coords[[1]])
        seqnames2 <- gsub(':.*', '', coords[[2]])
        df <- data.frame(
            chrom1 = seqnames1, 
            start1 = parsed_hic$x, 
            end1 = parsed_hic$x + resolution,
            chrom2 = seqnames2, 
            start2 = parsed_hic$y, 
            end2 = parsed_hic$y + resolution, 
            count = parsed_hic$counts, 
            score = parsed_hic_balanced$counts
        )
        df <- dplyr::filter(df, chrom1 == chrom2)
        df <- dplyr::mutate(df, 
            diag = {start2 - start1}/resolution, 
            distance = start2 - start1
        )
    }, error = function(e) {
        as.data.frame(InteractionSet::GInteractions())
    })}) |> dplyr::bind_rows()
    pixs$bin1_id <- left_join(
        pixs, bins, by = c(chrom1 = 'seqnames', end1 = 'end')
    )$bin_id
    pixs$bin2_id <- left_join(
        pixs, bins, by = c(chrom2 = 'seqnames', end2 = 'end')
    )$bin_id
    pixs <- dplyr::arrange(pixs, bin1_id, bin2_id) 
    pixs <- pixs[stats::complete.cases(pixs[ , c('bin1_id', 'bin2_id')]), ]
    res <- list(
        bins = bins, 
        pixels = pixs
    )

    return(res)
}
