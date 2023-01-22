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
            coords <- char2coords(coords)
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
            ean1 <- end(an1)
            maxean1 <- GenomeInfoDb::seqlengths(si)[GenomeInfoDb::seqlevels(si) == coords[[1]]]
            end(an1)[ean1 > maxean1] <- maxean1
            GenomeInfoDb::seqlevels(an1) <- GenomeInfoDb::seqlevels(si)
            GenomeInfoDb::seqinfo(an1) <- si
            an2 <- GenomicRanges::GRanges(
                seqnames = seqnames2, 
                ranges = IRanges::IRanges(
                    start = parsed_hic$y + 1, 
                    width = resolution
                )
            )
            ean2 <- end(an2)
            maxean2 <- GenomeInfoDb::seqlengths(si)[GenomeInfoDb::seqlevels(si) == coords[[2]]]
            end(an2)[ean2 > maxean2] <- maxean2
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
        an2 <- GenomicRanges::GRanges(
            seqnames = full_parsed_hic$seqnames2, 
            ranges = IRanges::IRanges(
                start = full_parsed_hic$start2, 
                width = resolution
            )
        )
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
        ean1 <- end(an1)
        maxean1 <- GenomeInfoDb::seqlengths(si)[GenomeInfoDb::seqlevels(si) == seqnames]
        end(an1)[ean1 > maxean1] <- maxean1
        GenomeInfoDb::seqlevels(an1) <- GenomeInfoDb::seqlevels(si)
        GenomeInfoDb::seqinfo(an1) <- si
        an2 <- GenomicRanges::GRanges(
            seqnames = seqnames, 
            ranges = IRanges::IRanges(
                start = parsed_hic$y + 1, 
                width = resolution
            )
        )
        ean2 <- end(an2)
        maxean2 <- GenomeInfoDb::seqlengths(si)[GenomeInfoDb::seqlevels(si) == seqnames]
        end(an2)[ean2 > maxean2] <- maxean2
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
        stop("Non-symetric parsing of .hic files is currently unsupported.")
    }

    # Match anchor ID with each interaction
    gi$bin_id1 <- S4Vectors::subjectHits(
        GenomicRanges::findOverlaps(an1, anchors)
    ) - 1
    gi$bin_id2 <- S4Vectors::subjectHits(
        GenomicRanges::findOverlaps(an2, anchors)
    ) - 1
    
    # Add extra info
    InteractionSet::regions(gi)$chr <- GenomicRanges::seqnames(InteractionSet::regions(gi))
    InteractionSet::regions(gi)$center <- GenomicRanges::start(GenomicRanges::resize(InteractionSet::regions(gi), fix = "center", width = 1))
    InteractionSet::regions(gi)$bin_id <- anchors$bin_id[BiocGenerics::match(regions(gi), anchors)]
    
    return(gi)
}

#' @param verbose Print resolutions in the console
#' @return vector
#'
#' @import rhdf5
#' @rdname parse-hic
#' @export
#' @examples 
#' hicPath <- HiContactsData::HiContactsData('yeast_wt', 'hic')
#' lsHicResolutions(hicPath, verbose = TRUE)

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
        dplyr::filter(seqnames != 'ALL') |>
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
