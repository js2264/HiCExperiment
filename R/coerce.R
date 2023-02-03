################################################################################
################################################################################
###############                                                  ###############
###############                    COERCING                      ###############
###############                                                  ###############
################################################################################
################################################################################

#' @name as
#' @export
#' @rdname HiCExperiment

setAs("HiCExperiment", "GInteractions", function(from) interactions(from))

#' @name as
#' @export
#' @rdname HiCExperiment

setAs("HiCExperiment", "InteractionSet", function(from) {
    gi <- interactions(from)
    is <- InteractionSet::InteractionSet(
        gi, 
        assays = lapply(
            names(scores(from)), function(n) as.matrix(scores(from, n))
        )
    )
    return(is)
})

#' @name as
#' @export
#' @rdname HiCExperiment

setAs("HiCExperiment", "ContactMatrix", function(from) {
    x <- interactions(from, fillout.regions = TRUE)
    if ('balanced' %in% names(scores(from))) {
        x$score <- scores(from, 'balanced')
        gi2cm(x)
    } 
    else {
        x$score <- scores(from, 1)
        gi2cm(x)
    }
})

#' @name as
#' @export
#' @rdname HiCExperiment

setAs("HiCExperiment", "matrix", function(from) {
    as(from, "ContactMatrix") |> cm2matrix(sparse = TRUE)
})

#' @name as
#' @export
#' @rdname HiCExperiment

setAs("HiCExperiment", "data.frame", function(from) {
    x <- interactions(from)
    x <- as.data.frame(x)
    x <- x[, !colnames(x) %in% c("chr1", "chr2", "bin_id1.1", "bin_id2.1")]
    for (n in names(scores(from))) {
        x[[n]] <- scores(from, n)
    }
    return(x)
})

#' @name as
#' @export
#' @rdname HiCExperiment

setMethod("as.matrix", "HiCExperiment", function(x) {
    xx <- as(x, 'matrix')
})

#' @name as
#' @export
#' @rdname HiCExperiment

setMethod("as.data.frame", "HiCExperiment", function(x) {
    as(x, 'data.frame')
})


################################################################################
################################################################################
###############                                                  ###############
###############              OTHER COERCING                      ###############
###############                                                  ###############
################################################################################
################################################################################

#' Other coercing functions
#' 
#' @import InteractionSet
#' @importFrom GenomicRanges mcols
#' @param gi GInteractions object
#' @param use.scores Which scores to use to inflate GInteractions
#' @rdname parse-other
#' @export
#' @examples 
#' mcoolPath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' contacts <- import(mcoolPath, focus = 'XVI', resolution = 16000, format = 'cool')
#' gis <- interactions(contacts)
#' cm <- gi2cm(gis, 'balanced')
#' cm

gi2cm <- function(gi, use.scores = 'score') {
    if (!use.scores %in% colnames(S4Vectors::mcols(gi))) 
        stop("`use.scores` argument not found in the provided interactions")
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
#' @param seqnames1,start1,end1,seqnames2,start2,end2 Names (as strings) of 
#' columns containing corresponding information in a data.frame parsed into 
#' GInteractions
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
    if (!is.na(replace_NA)) m[is.na(m)] <- replace_NA
    if (!sparse) m <- base::as.matrix(m)
    m
}

#' @param df A `data.frame` object
#' @return a `GInteractions` object
#'
#' @rdname parse-other
#' @export
#' @examples 
#' df2gi(data.frame(
#'     chr1 = 'I', start1 = 10, end1 = 100, 
#'     chr2 = 'I', start2 = 40, end2 = 1000, 
#'     score = 12, 
#'     weight = 0.234, 
#'     filtered = TRUE
#' ), seqnames1 = 'chr1', seqnames2 = 'chr2')

df2gi <- function(df, seqnames1 = 'seqnames1', start1 = 'start1', end1 = 'end1', seqnames2 = 'seqnames2', start2 = 'start2', end2 = 'end2') {
    gi <- InteractionSet::GInteractions(
        anchor1 = GenomicRanges::GRanges(
            df[[seqnames1]], IRanges::IRanges(df[[start1]], df[[end1]])
        ),
        anchor2 = GenomicRanges::GRanges(
            df[[seqnames2]], IRanges::IRanges(df[[start2]], df[[end2]])
        )
    )
    S4Vectors::mcols(gi) <- dplyr::select(
        df, !dplyr::any_of(c(seqnames1, start1, end1, seqnames2, start2, end2))
    )
    return(gi)
}

