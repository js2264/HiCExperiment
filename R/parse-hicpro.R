#' Parsing hicpro files (matrix & bed)
#' 
#' These functions are the workhorse internal functions used to import 
#' HiC-Pro `.matrix` and `.bed` files as GenomicInteractions (wrapped into a `HiCExperiment` object
#' by `HiCExperiment()` function).
#'
#' @param file file
#' @param bed bed
#' @return a GenomicInteractions object
#'
#' @import InteractionSet
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges resize
#' @rdname parse-hicpro

.hicpro2gi <- function(file, bed) {
    
    file <- gsub('~', Sys.getenv('HOME'), file)
    
    # Get raw counts for bins from hic
    matrix_df <- vroom::vroom(
        file, 
        col_names = FALSE, 
        progress = FALSE, 
        show_col_types = FALSE
    )
    colnames(matrix_df) <- c("start_idx", "stop_idx", "value")
    
    # Get anchors from hicpro
    anchors <- .getHicproAnchors(bed)
    an1 <- left_join(
        matrix_df, as.data.frame(anchors), by = c('start_idx' = 'bin_id')
    ) |> as("GRanges")
    GenomicRanges::mcols(an1) <- NULL
    an2 <- left_join(
        matrix_df, as.data.frame(anchors), by = c('stop_idx' = 'bin_id')
    ) |> as("GRanges")
    GenomicRanges::mcols(an2) <- NULL
    gi <- InteractionSet::GInteractions(
        an1, 
        an2, 
    )
    GenomeInfoDb::seqlevels(gi) <- GenomeInfoDb::seqlevels(anchors)
    GenomeInfoDb::seqinfo(gi) <- GenomeInfoDb::seqinfo(anchors)

    # Find bin IDs
    gi$bin_id1 <- S4Vectors::subjectHits(
        GenomicRanges::findOverlaps(an1, anchors)
    ) - 1
    gi$bin_id2 <- S4Vectors::subjectHits(
        GenomicRanges::findOverlaps(an2, anchors)
    ) - 1

    # Associate counts for bins to corresponding anchors
    gi$count <- matrix_df[[3]]
    
    # Add extra info
    InteractionSet::regions(gi)$chr <- GenomicRanges::seqnames(InteractionSet::regions(gi))
    InteractionSet::regions(gi)$center <- GenomicRanges::start(GenomicRanges::resize(InteractionSet::regions(gi), fix = "center", width = 1))
    InteractionSet::regions(gi)$bin_id <- anchors$bin_id[BiocGenerics::match(regions(gi), anchors)]
    
    return(gi)
}

#' @rdname parse-hicpro

.getHicproAnchors <- function(bed) {
    si <- .hicpro2seqinfo(bed)
    bed1 <- vroom::vroom(
        bed, 
        col_names = FALSE, 
        progress = FALSE, 
        show_col_types = FALSE, 
        n_max = 10
    )
    resolution <- max(unique(bed1[[3]][1] - bed1[[2]][1]))
    anchors <- GenomicRanges::tileGenome(
        si, tilewidth = resolution, cut.last.tile.in.chrom = TRUE
    )
    anchors$bin_id <- seq_along(anchors) - 1
    names(anchors) <- paste(GenomicRanges::seqnames(anchors), GenomicRanges::start(anchors), GenomicRanges::end(anchors), sep = "_")
    return(anchors)
}

#' @rdname parse-hicpro

.hicpro2seqinfo <- function(bed) {
    anchors_df <- vroom::vroom(
        bed, 
        col_names = FALSE, 
        progress = FALSE, 
        show_col_types = FALSE
    )
    anchors_df |> 
        group_by(X1) |> 
        summarize(max = max(X3)) |> 
        dplyr::rename(seqnames = X1, seqlengths = max) |> 
        as.data.frame() |> 
        as("Seqinfo")
}