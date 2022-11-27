#' Parsing hic files
#' 
#' These functions are the workhorse internal functions used to import 
#' a `.hic` file as GenomicInteractions (wrapped into a `HiCExperiment` object
#' by `HiCExperiment()` function).
#'
#' @param file file
#' @param coords NULL, character, or GRanges. 
#'   Can also be a Pairs object of paired GRanges (length of 1).
#' @param resolution resolution
#' @return a GenomicInteractions object
#'
#' @import InteractionSet
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges resize
#' @rdname parse-hic

.hic2gi <- function(file, coords = NULL, resolution = NULL) {
   
    file <- gsub('~', Sys.getenv('HOME'), file)
    
    # Mutate Pairs provided as characters to real Pairs
    if (!is.null(coords)) {
        if (grepl(' x ', coords)) {
            coords <- char2coords(coords)
        }
    }

    # Check if the provided coords are GRanges or Pairs
    is_pair <- is(coords, 'Pairs')

    # Get raw counts for bins from hic
    if (!is_pair) {
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
    }
    else {
        stop("unsupported xx")
    }
    an1 <- GenomicRanges::GRanges(
        seqnames = seqnames, 
        ranges = IRanges::IRanges(
            start = parsed_hic$x + 1, 
            width = resolution
        )
    )
    an2 <- GenomicRanges::GRanges(
        seqnames = seqnames, 
        ranges = IRanges::IRanges(
            start = parsed_hic$y + 1, 
            width = resolution
        )
    )
    gi <- InteractionSet::GInteractions(
        an1, 
        an2, 
    )

    # Find bin IDs
    si <- strawr::readHicChroms(file) |> 
        setNames(c('seqnames', 'seqlengths')) |>
        as("Seqinfo")
    anchors <- GenomicRanges::tileGenome(
        si, tilewidth = resolution, cut.last.tile.in.chrom = TRUE
    )
    anchors$bin_id <- seq_along(anchors)
    gi$bin_id1 <- S4Vectors::subjectHits(
        GenomicRanges::findOverlaps(an1, anchors)
    )
    gi$bin_id2 <- S4Vectors::subjectHits(
        GenomicRanges::findOverlaps(an2, anchors)
    )

    # Associate raw counts for bins to corresponding anchors
    gi$count <- parsed_hic$counts
    gi$score <- parsed_hic_balanced$counts
    
    # Add extra info
    InteractionSet::regions(gi)$chr <- GenomicRanges::seqnames(InteractionSet::regions(gi))
    InteractionSet::regions(gi)$center <- GenomicRanges::start(GenomicRanges::resize(InteractionSet::regions(gi), fix = "center", width = 1))
    InteractionSet::regions(gi)$bin_id <- anchors$bin_id[BiocGenerics::match(regions(gi), anchors)]
    
    return(gi)
}

#' @param file file
#' @param verbose verbose
#' @return vector
#'
#' @import rhdf5
#' @import stringr
#' @import tidyr
#' @import dplyr
#' @import GenomicInteractions
#' @rdname parse

.lsHicFiles <- function(file, verbose = FALSE) {
    tryCatch(
        expr = {strawr::readHicChroms(file); TRUE}, 
        error = function(e) {message("Provided file is not a .hic file.")}
    )
}

#' @param file file
#' @param verbose Print resolutions in the console
#' @return vector
#'
#' @import rhdf5
#' @import tidyr
#' @rdname parse
#' @export

lsHicResolutions <- function(file, verbose = FALSE) {
    res <- strawr::readHicBpResolutions(file)
    if (verbose) message(S4Vectors::coolcat("resolutions(%d): %s", res))
    invisible(as.integer(res))
}

