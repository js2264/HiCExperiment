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
#' contacts <- import(mcool_path, format = 'cool')
#' gis <- interactions(contacts)
#' gis$score <- scores(contacts, 1)
#' cm <- gi2cm(gis)
#' cm

gi2cm <- function(gi, use.scores = 'score') {
    InteractionSet::inflate(
        gi,
        rows = seq_along(InteractionSet::regions(gi)),
        columns = seq_along(InteractionSet::regions(gi)),
        fill = GenomicRanges::mcols(gi)[[use.scores]]
    )
}

#' @param cm A `ContactMatrix` object
#' @param replace_NA Replace NA values
#' @return a dense matrix
#'
#' @importFrom Matrix as.matrix
#' @rdname parse-other
#' @export
#' @examples 
#' cm2matrix(cm)[1:10, 1:10]

cm2matrix <- function(cm, replace_NA = NA) {
    m <- Matrix::as.matrix(cm)
    m[is.na(m)] <- replace_NA
    m
}
