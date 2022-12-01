#' Other parsing functions
#' 
#' @import InteractionSet
#' @importFrom GenomicRanges mcols
#' @param gi GInteractions
#' @rdname parse-other
#' @export

gi2cm <- function(gi) {
    InteractionSet::inflate(
        gi,
        rows = seq_along(InteractionSet::regions(gi)),
        columns = seq_along(InteractionSet::regions(gi)),
        fill = GenomicRanges::mcols(gi)[['score']]
    )
}

#' @param cm A `ContactMatrix` object
#' @param replace_NA Replace NA values
#' @return a dense matrix
#'
#' @importFrom Matrix as.matrix
#' @rdname parse-other
#' @export

cm2matrix <- function(cm, replace_NA = NA) {
    m <- Matrix::as.matrix(cm)
    m[is.na(m)] <- replace_NA
    m
}
