#' Other parsing functions
#' 
#' @import InteractionSet
#' @importFrom GenomicRanges mcols
#' @param gi GInteractions
#' @param use.scores Which scores to use to inflate GInteractions
#' @rdname parse-other
#' @export
#' @examples 
#' mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' contacts <- import(mcool_path, focus = 'XVI', resolution = 16000, format = 'cool')
#' gis <- interactions(contacts)
#' cm <- gi2cm(gis, 'balanced')
#' cm

gi2cm <- function(gi, use.scores = 'score') {
    # if (!use.scores %in% colnames(S4Vectors::mcols(gi))) 
    #     stop("`use.scores` argument not found in the provided interactions")
    # seqnames <- unique(GenomeInfoDb::seqnames(regions(gi)))
    # if (length(seqnames == 1)) {
    #     gr <- GenomicRanges::GRanges(seqinfo(gi))
    #     res <- GenomicRanges::width(regions(gi))[1]
    #     gr <- GenomicRanges::tileGenome(
    #         seqinfo(gi), 
    #         tilewidth = res, 
    #         cut.last.tile.in.chrom = TRUE
    #     )
    #     gr <- gr[GenomicRanges::seqnames(gr) == seqnames]
    #     InteractionSet::replaceRegions(gi) <- gr
    # }
    an <- anchors(gi)
    allTrans <- all(GenomicRanges::seqnames(an[[1]]) != GenomicRanges::seqnames(an[[2]]))
    if (allTrans) {
        stop("Only trans contacts are found.")
    }
    InteractionSet::inflate(
        gi,
        rows = InteractionSet::regions(gi), 
        columns = InteractionSet::regions(gi), 
        fill = GenomicRanges::mcols(gi)[[use.scores]], 
        sparse = TRUE
    )
}

#' @param cm A `ContactMatrix` object
#' @param replace_NA Replace NA values
#' @param sparse Whether to return the contact matrix as a sparse matrix 
#' (default: FALSE)
#' @return a dense matrix
#'
#' @importFrom Matrix as.matrix
#' @rdname parse-other
#' @export
#' @examples 
#' cm2matrix(cm)[1:10, 1:10]

cm2matrix <- function(cm, replace_NA = NA, sparse = FALSE) {
    m <- Matrix::as.matrix(cm)
    m <- as(m, 'TsparseMatrix')
    m[is.na(m)] <- replace_NA
    if (!sparse) m <- base::as.matrix(m)
    m
}
