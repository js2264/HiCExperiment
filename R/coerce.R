#' @title Coercing functions
#' @name as
#' @rdname as
#' 
#' @description
#' Coercing functions available for HiCExperiment objects.
#' 
#' @aliases coerce,HiCExperiment,GInteractions-method
#' @aliases coerce,HiCExperiment,InteractionSet-method
#' @aliases coerce,HiCExperiment,ContactMatrix-method
#' @aliases coerce,HiCExperiment,matrix-method
#' @aliases coerce,HiCExperiment,data.frame-method
#' @aliases as.matrix,HiCExperiment-method
#' @aliases as.data.frame,HiCExperiment-method
#' 
#' @import InteractionSet
#' @importFrom GenomicRanges mcols
#' @param gi GInteractions object
#' @param x HiCExperiment object
#' @param use.scores Which scores to use to inflate GInteractions
#' @param cm A `ContactMatrix` object
#' @param replace_NA Replace NA values
#' @param sparse Whether to return the contact matrix as a sparse matrix 
#' @param seqnames1,start1,end1,seqnames2,start2,end2 Names (as strings) of 
#' columns containing corresponding information in a data.frame parsed into 
#' GInteractions
#' (default: FALSE)
#' @param df A `data.frame` object
#' @importFrom Matrix as.matrix
#' @export
#' @examples 
#' mcoolPath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' contacts <- import(mcoolPath, focus = 'XVI', resolution = 16000, format = 'cool')
#' gis <- interactions(contacts)
#' cm <- gi2cm(gis, 'balanced')
#' cm
#' cm2matrix(cm)[1:10, 1:10]
#' df2gi(data.frame(
#'     chr1 = 'I', start1 = 10, end1 = 100, 
#'     chr2 = 'I', start2 = 40, end2 = 1000, 
#'     score = 12, 
#'     weight = 0.234, 
#'     filtered = TRUE
#' ), seqnames1 = 'chr1', seqnames2 = 'chr2')
NULL 

################################################################################
################################################################################
###############                                                  ###############
###############     HICEXPERIMENT COERCING                       ###############
###############                                                  ###############
################################################################################
################################################################################

#' @export
#' @name as

setAs("HiCExperiment", "GInteractions", function(from) interactions(from))

#' @export
#' @name as

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

#' @export
#' @name as

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

#' @export
#' @name as

setAs("HiCExperiment", "matrix", function(from) {
    as(from, "ContactMatrix") |> cm2matrix(sparse = TRUE)
})

#' @export
#' @name as

setAs("HiCExperiment", "data.frame", function(from) {
    x <- interactions(from)
    x <- as.data.frame(x)
    x <- x[, !colnames(x) %in% c("chr1", "chr2", "bin_id1.1", "bin_id2.1")]
    for (n in names(scores(from))) {
        x[[n]] <- scores(from, n)
    }
    return(x)
})

#' @export
#' @name as

setMethod("as.matrix", "HiCExperiment", function(x) {
    xx <- as(x, 'matrix')
})

#' @export
#' @name as

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

#' @export
#' @name as

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

#' @name as
#' @export

cm2matrix <- function(cm, replace_NA = NA, sparse = FALSE) {
    m <- Matrix::as.matrix(cm)
    m <- as(m, 'TsparseMatrix')
    if (!is.na(replace_NA)) m[is.na(m)] <- replace_NA
    if (!sparse) m <- base::as.matrix(m)
    m
}

#' @name as
#' @export

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

