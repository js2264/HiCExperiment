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

#' @param file pairs file: `<readname>\t<chr1>\t<start1>\t<chr2>\t<start2>`
#' @param chr1.field chr1.field
#' @param start1.field start1.field
#' @param chr2.field chr2.field
#' @param start2.field start2.field
#' @param strand1.field strand1.field
#' @param strand2.field strand2.field
#' @param frag1.field frag1.field
#' @param frag2.field frag2.field
#' @param nThread Number of CPUs to use to import the `pairs` file in R
#' @param nrows Number of pairs to import
#' @return a GenomicInteractions object
#'
#' @importFrom vroom vroom
#' @importFrom glue glue
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicInteractions GenomicInteractions
#' @importFrom GenomicInteractions calculateDistances
#' @importFrom IRanges IRanges
#' @import tibble
#' @rdname parse-other

pairs2gi <- function(
    file, 
    chr1.field = 2, 
    start1.field = 3, 
    chr2.field = 4, 
    start2.field = 5, 
    strand1.field = 6, 
    strand2.field = 7, 
    frag1.field = NULL, 
    frag2.field = NULL, 
    nThread = 16, 
    nrows = Inf  
) {
    sel1 <- dplyr::all_of(c(chr1.field, start1.field, strand1.field, frag1.field))
    sel2 <- dplyr::all_of(c(chr2.field, start2.field, strand2.field, frag2.field))
    anchors1 <- vroom::vroom(
        file,
        n_max = nrows,
        col_select = sel1,
        comment = '#',
        col_names = FALSE,
        show_col_types = FALSE
    )
    anchors2 <- vroom::vroom(
        file,
        n_max = nrows,
        col_select = sel2,
        comment = '#',
        col_names = FALSE,
        show_col_types = FALSE
    )  
    anchor_one <- GenomicRanges::GRanges(
        anchors1[[1]],
        IRanges::IRanges(anchors1[[2]], width = 1), 
        strand = anchors1[[3]]
    )
    anchor_two <- GenomicRanges::GRanges(
        anchors2[[1]],
        IRanges::IRanges(anchors2[[2]], width = 1), 
        strand = anchors2[[3]]
    )
    gi <- GenomicInteractions::GenomicInteractions(anchor_one, anchor_two)
    if (!is.null(frag1.field) & !is.null(frag2.field)) {
        gi$frag1 <- anchors1[[4]]
        gi$frag2 <- anchors2[[4]]
    } 
    else {
        gi$frag1 <- NA
        gi$frag2 <- NA
    }
    gi$distance <- GenomicInteractions::calculateDistances(gi) 
    return(gi)
}
