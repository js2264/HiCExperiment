#' @title HiCExperiment binning methods
#' 
#' @name bin-methods
#' @aliases bin,GInteractions,character-method
NULL

#' @rdname bin-methods
#' @export 
#' @examples 
#' pairsf <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#' pf <- PairsFile(pairsf)
#' 

setMethod("bin", signature(x = "GInteractions"), function(x, resolution, seqinfo = NULL, pairsFile = NULL) {
    
    ## -- Approximate seqinfo if missing
    if (is.null(seqinfo)) {
        # sn <- unique(c(levels(seqnames(an[[1]]))), c(levels(seqnames(an[[2]]))))
        # sl <- NULL
        # seqinfo <- GenomeInfoDb::Seqinfo(
        #     seqnames = unique(c())
        # )
        stop("Currently not supported")
    }

    ## -- Define regions
    re <- GenomicRanges::tileGenome(
        seqinfo, tilewidth = resolution, cut.last.tile.in.chrom = TRUE
    )
    re$bin_id <- seq_along(re) - 1

    ## -- create ubinnedgi
    an <- anchors(x)
    print("HI")
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
    hic <- HiContacts::normalize(hic)
    names(hic@scores)[2] <- 'balanced'
    methods::validObject(hic)
    return(hic)
})


setMethod("bin", signature(x = "PairsFile"), function(x, resolution, seqinfo = NULL) {
    
    ## -- Import pairs as GInteractions
    gi <- import(x)
    gi$frag1 <- NULL
    gi$frag2 <- NULL
    gi$distance <- NULL

    bin(gi, resolution, seqinfo, pairsFile = path(x))

})
