#' @title HiCExperiment binning methods
#' 
#' @name bin-methods
#' @aliases bin,GInteractions,character-method
#' @param x A PairsFile or GInteractions object
#' @param resolution Which resolution to use to bin the interactions
#' @param seqinfo Seqinfo object
#' @examples 
#' pairsf <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#' pf <- PairsFile(pairsf)
NULL

#' @rdname bin-methods

setMethod("bin", signature(x = "GInteractions", resolution = "numeric"), function(x, resolution, seqinfo = NULL) {

    ## -- Approximate seqinfo if missing
    if (is.null(seqinfo)) {
        re_df <- regions(gi) |> 
            as.data.frame() |> 
            dplyr::group_by(seqnames) |> 
            dplyr::summarize(seqlengths = max(end))
        seqinfo <- GenomeInfoDb::Seqinfo(
            seqnames = levels(re_df$seqnames), 
            seqlengths = re_df$seqlengths
        )
    }

    ## -- Define regions
    re <- GenomicRanges::tileGenome(
        seqinfo, tilewidth = resolution, cut.last.tile.in.chrom = TRUE
    )
    re$bin_id <- seq_along(re) - 1

    ## -- create ubinnedgi
    an <- anchors(x)
    binnedgi <- InteractionSet::GInteractions(
        S4Vectors::subjectHits(InteractionSet::findOverlaps(an[[1]], re)), 
        S4Vectors::subjectHits(InteractionSet::findOverlaps(an[[2]], re)), 
        regions = re
    )
    ubinnedgi <- unique(binnedgi)
    uan <- anchors(ubinnedgi)
    ubinnedgi$bin_id1 <- uan[[1]]$bin_id
    ubinnedgi$bin_id2 <- uan[[2]]$bin_id
    cnts <- as.data.frame(binnedgi) |> 
        dplyr::group_by(seqnames1, start1, seqnames2, start2) |> 
        dplyr::tally() |> 
        dplyr::pull(n)

    ## -- create binned HiCExperiment object
    pairsFile <- attr(gi, 'pairs')[[1]]
    hic <- methods::new("HiCExperiment", 
        fileName = ifelse(!is.null(pairsFile), pairsFile, ""),
        focus = NULL, 
        resolutions = resolution, 
        resolution = resolution, 
        interactions = ubinnedgi, 
        scores = S4Vectors::SimpleList(
            'count' = as.numeric(cnts)
        ), 
        topologicalFeatures = S4Vectors::SimpleList(), 
        pairsFile = ifelse(!is.null(pairsFile), pairsFile, ""), 
        metadata = list(
            note = "Binned pairs file"
        )
    )

    ## -- Balance HiCExperiment 
    if (requireNamespace("HiContacts", quietly = TRUE)) {
        hic <- HiContacts::normalize(hic)
    }
    else {
        warning('Install `HiContacts` package (`BiocManager::install("HiContacts")`)\nto balance Hi-C data.')
    }
    names(hic@scores)[2] <- 'balanced'

    ## -- Return new HiCExperiment
    methods::validObject(hic)
    return(hic)
})

#' @rdname bin-methods

setMethod("bin", signature(x = "PairsFile", resolution = "numeric"), function(x, resolution, seqinfo = NULL) {
    ## -- Import pairs as GInteractions
    gi <- import(x)
    gi$frag1 <- NULL
    gi$frag2 <- NULL
    gi$distance <- NULL
    attr(gi, 'pairs') <- path(x)
    bin(gi, resolution, seqinfo)
})
