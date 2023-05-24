#' @title `HiCExperiment` methods
#' 
#' @name HiCExperiment
#' @rdname HiCExperiment
#' @aliases resolutions,HiCExperiment-method
#' @aliases resolution,HiCExperiment-method
#' @aliases focus,HiCExperiment-method
#' @aliases focus<-,HiCExperiment-method
#' @aliases focus<-,HiCExperiment,character-method
#' @aliases zoom,HiCExperiment,numeric-method
#' @aliases refocus,HiCExperiment,character-method
#' @aliases scores,HiCExperiment-method
#' @aliases scores,HiCExperiment,missing-method
#' @aliases scores,HiCExperiment,character-method
#' @aliases scores,HiCExperiment,numeric-method
#' @aliases scores<-,HiCExperiment-method
#' @aliases scores<-,HiCExperiment,character,numeric-method
#' @aliases topologicalFeatures,HiCExperiment-method
#' @aliases topologicalFeatures,HiCExperiment,missing-method
#' @aliases topologicalFeatures,HiCExperiment,character-method
#' @aliases topologicalFeatures,HiCExperiment,numeric-method
#' @aliases topologicalFeatures<-,HiCExperiment-method
#' @aliases topologicalFeatures<-,HiCExperiment,character,GRangesOrGInteractions-method
#' @aliases pairsFile,HiCExperiment-method
#' @aliases pairsFile<-,HiCExperiment-method
#' @aliases pairsFile<-,HiCExperiment,character-method
#' @aliases metadata<-,HiCExperiment-method
#' @aliases metadata<-,HiCExperiment,list-method
#' @aliases fileName,HiCExperiment-method
#' @aliases interactions,HiCExperiment-method
#' @aliases interactions<-,HiCExperiment-method
#' @aliases length,HiCExperiment-method
#' @aliases [,HiCExperiment-method
#' @aliases [,HiCExperiment,numeric,ANY,ANY-method
#' @aliases subsetByOverlaps,HiCExperiment,GRanges-method
#' @aliases [,HiCExperiment,logical,ANY,ANY-method
#' @aliases [,HiCExperiment,character,ANY,ANY-method
#' @aliases [,HiCExperiment,GInteractions,ANY,ANY-method
#' @aliases [,HiCExperiment,Pairs,ANY,ANY-method
#' @aliases seqinfo,HiCExperiment-method
#' @aliases bins,HiCExperiment-method
#' @aliases anchors,HiCExperiment-method
#' @aliases regions,HiCExperiment-method
#' @aliases cis,HiCExperiment-method
#' @aliases trans,HiCExperiment-method
#' @aliases show,HiCExperiment-method
#' 
#' @param x A \code{HiCExperiment} object.
#' @param object A \code{HiCExperiment} object.
#' @param name Name of the element to access in topologicalFeatures or scores SimpleLists.
#' @param value Value to add to topologicalFeatures, scores, pairsFile or metadata slots.
#' @param i,ranges a GRanges, coordinates in character, or boolean vector to subset a HiCExperiment
#' @param type any of `within` or `any`, to subset interactions by overlap with a provided GRanges.
#' @param fillout.regions Whehter to add missing regions to GInteractions' regions? 
#' 
#' @importMethodsFrom BiocGenerics fileName
#' @importFrom GenomeInfoDb seqinfo
#' 
#' @include AllGenerics.R
NULL

################################################################################
################################################################################
###############                                                  ###############
###############          METHODS FOR NEW GENERICS                ###############
###############                                                  ###############
################################################################################
################################################################################

#' @export
#' @rdname HiCExperiment

setMethod("resolutions", "HiCExperiment", function(x) x@resolutions)

#' @export
#' @rdname HiCExperiment

setMethod("resolution", "HiCExperiment", function(x) x@resolution)

#' @export
#' @rdname HiCExperiment

setMethod("focus", "HiCExperiment", function(x) x@focus)

#' @export
#' @rdname HiCExperiment

setMethod("focus<-", signature(x = "HiCExperiment", value = "character"), function(x, value) {
    x@focus <- value
    x
})

#' @export
#' @rdname HiCExperiment

setMethod("zoom", c("HiCExperiment", "numeric"), function(x, resolution) {
    HiCExperiment(
        fileName(x), 
        resolution = as.integer(resolution), 
        focus = focus(x), 
        metadata = S4Vectors::metadata(x), 
        topologicalFeatures = topologicalFeatures(x), 
        pairsFile = pairsFile(x)
    )
})

#' @export
#' @rdname HiCExperiment

setMethod("refocus", c("HiCExperiment", "character"), function(x, focus) {
    if (.is_cool(fileName(x))) {
        res <- NULL
    } else {
        res <- resolution(x)
    }
    HiCExperiment(
        fileName(x), 
        resolution = res, 
        focus = focus, 
        metadata = S4Vectors::metadata(x), 
        topologicalFeatures = topologicalFeatures(x), 
        pairsFile = pairsFile(x)
    )
})

#' @export
#' @rdname HiCExperiment

setMethod("scores", signature(x = "HiCExperiment", name = "missing"), function(x) x@scores)

#' @export
#' @rdname HiCExperiment

setMethod("scores", signature(x = "HiCExperiment", name = "character"), function(x, name) {
    .check_scores(x, name)
    return(x@scores[[name]])
})

#' @export
#' @rdname HiCExperiment

setMethod("scores", signature(x = "HiCExperiment", name = "numeric"), function(x, name) {
    if (name > length(scores(x))) {
        stop(paste0('Only ', length(scores(x)), ' scores in x.'))
    }
    return(x@scores[[name]])
})

#' @export
#' @rdname HiCExperiment

setMethod("scores<-", c(x = "HiCExperiment", name = "character", value = "numeric"), function(x, name, value) {
    x@scores[[name]] <- value
    return(x)
})

#' @export
#' @rdname HiCExperiment

setMethod("topologicalFeatures", signature(x = "HiCExperiment", name = "missing"), function(x) {
    S4Vectors::SimpleList(as.list(x@topologicalFeatures))
})

#' @export
#' @rdname HiCExperiment

setMethod("topologicalFeatures", signature(x = "HiCExperiment", name = "character"), function(x, name) {
    if (!name %in% names(topologicalFeatures(x))) {
        stop(paste0(name, ' not in topologicalFeatures.'))
    }
    x@topologicalFeatures[[name]]
})

#' @export
#' @rdname HiCExperiment

setMethod("topologicalFeatures", signature(x = "HiCExperiment", name = "numeric"), function(x, name) {
    if (name > length(topologicalFeatures(x))) {
        stop(paste0('Only ', length(topologicalFeatures(x)), ' topologicalFeatures in x.'))
    }
    x@topologicalFeatures[[name]]
})

#' @export
#' @rdname HiCExperiment

setMethod("topologicalFeatures<-", signature(x = "HiCExperiment", name = "character", value = "GRangesOrGInteractions"), function(x, name, value) {
    x@topologicalFeatures[[name]] <- value
    return(x)
})

#' @export
#' @rdname HiCExperiment

setMethod("pairsFile", "HiCExperiment", function(x) {
    x@pairsFile
})

#' @export
#' @rdname HiCExperiment

setMethod("pairsFile<-", signature(x = "HiCExperiment", value = "character"), function(x, value) {
    if (!file.exists(value)) {
        stop("Provided pairsFile does not exist. Aborting now.")
    }
    x@pairsFile <- value
    x
})

#' @export
#' @rdname HiCExperiment

setMethod("metadata<-", signature(x = "HiCExperiment", value = "list"), function(x, value) {
    x@metadata <- value
    x
})

################################################################################
################################################################################
###############                                                  ###############
###############               METHODS FOR SUBSET                 ###############
###############                                                  ###############
################################################################################
################################################################################

#' @export
#' @rdname HiCExperiment

setMethod("subsetByOverlaps", signature = c(x="HiCExperiment", ranges="numeric"), function(x, ranges) {
    subints <- interactions(x)[ranges]
    subre <- IRanges::subsetByOverlaps(regions(x), subints)
    InteractionSet::replaceRegions(subints) <- subre
    interactions(x) <- subints
    for (n in names(scores(x))) {
        scores(x, n) <- scores(x, n)[ranges]
    }
    return(x)
})

#' @export
#' @rdname HiCExperiment

setMethod("subsetByOverlaps", signature = c(x="HiCExperiment", ranges="logical"), function(x, ranges) {
    # Redirect to subsetByOverlaps,HiCExperiment,numeric-method
    x[which(ranges)]
})

#' @export
#' @rdname HiCExperiment

setMethod("subsetByOverlaps", signature = c(x="HiCExperiment", ranges="GRanges"), function(x, ranges, type = c('within', 'any')) {
    type <- match.arg(type)
    if (type == "any") { ## Anchors not constrained within segment
        i <- IRanges::overlapsAny(
            interactions(x), ranges, type = 'within'
        )
        # Redirect to subsetByOverlaps,HiCExperiment,numeric-method
        return(x[i])
    }
    else if (type == "within") { ## Anchors constrained within segment, for each range
        # Redirect to subsetByOverlaps,HiCExperiment,GInteractions-method
        subsetByOverlaps(x, InteractionSet::GInteractions(ranges, ranges))
    }
    else {
        stop("'type' should be one of \"any\", \"both\"")
    }
})

#' @export
#' @rdname HiCExperiment

setMethod("subsetByOverlaps", signature = c(x="HiCExperiment", ranges="GInteractions"), function(x, ranges) {
    ranges <- swapAnchors(ranges)
    ints <- interactions(x)
    an <- anchors(ints)
    i <- lapply(seq_along(ranges), function(k) {
        which(IRanges::overlapsAny(
            an[["first"]], anchors(ranges[k], "first"), type = 'within'
        ) & IRanges::overlapsAny(
            an[["second"]], anchors(ranges[k], "second"), type = 'within'
        ))
    }) |> unlist() |> unique()
    # Redirect to subsetByOverlaps,HiCExperiment,numeric-method
    return(x[i])
})

#' @export
#' @rdname HiCExperiment

setMethod("subsetByOverlaps", signature = c(x="HiCExperiment", ranges="Pairs"), function(x, ranges) {
    gi <- GInteractions(S4Vectors::first(ranges), S4Vectors::second(ranges))
    # Redirect to subsetByOverlaps,HiCExperiment,GInteractions-method
    return(x[gi])
})

#' @export
#' @rdname HiCExperiment

setMethod("[", signature("HiCExperiment", "numeric"), function(x, i) {
    # Redirect to subsetByOverlaps,HiCExperiment,numeric-method
    subsetByOverlaps(x, i)
})

#' @export
#' @rdname HiCExperiment

setMethod("[", signature("HiCExperiment", "GRanges"), function(x, i) {
    # Redirect to subsetByOverlaps,HiCExperiment,GRanges-method
    subsetByOverlaps(x, i)
})

#' @export
#' @rdname HiCExperiment

setMethod("[", signature("HiCExperiment", "logical"), function(x, i) {
    # Redirect to subsetByOverlaps,HiCExperiment,logical-method
    subsetByOverlaps(x, i)
})

#' @export
#' @rdname HiCExperiment

setMethod("[", signature("HiCExperiment", "GInteractions"), function(x, i) {
    # Redirect to subsetByOverlaps,HiCExperiment,GInteractions-method
    subsetByOverlaps(x, i)
})

#' @export
#' @rdname HiCExperiment

setMethod("[", signature("HiCExperiment", "Pairs"), function(x, i) {
    # Redirect to subsetByOverlaps,HiCExperiment,Pairs-method
    subsetByOverlaps(x, i)
})

#' @export
#' @rdname HiCExperiment

setMethod("[", signature("HiCExperiment", "character"), function(x, i) {

    re_ <- regions(x)
    ints_ <- interactions(x)

    # Different possible situations: 
    #   'II:10001-20000'                    --> Redirects to GRanges method
    #   'II'                                --> Redirects to GRanges method
    #   'II:10001-20000|III:50001-90000'    --> Redirects to GInteractions method
    #   'II|III'                            --> Redirects to GInteractions method
    #   c('II', 'III')                      --> Redirects to GRanges method

    if (length(i) == 1) { # 'II:10001-20000', 'II:10001-20000|III:50000-90000', 'II' or 'II|III'
        if (
            grepl(
                '.*:[0-9]*-[0-9]*\\|[A-Za-z0-9]*:[0-9]*-[0-9]*$', i
            ) | grepl(
                '.*:[0-9]*-[0-9]*$', i
            ) 
        ) { # e.g. 'II:10001-20000' or 'II:10001-20000|III:50000-90000'
            i_ <- char2coords(i)
            valid_regions_first <- subsetByOverlaps(
                re_, S4Vectors::first(i_), type = 'within'
            )$bin_id
            valid_regions_second <- subsetByOverlaps(
                re_, S4Vectors::second(i_), type = 'within'
            )$bin_id
        }
        else if (
            grepl(
                '[^:]*\\|[^:]*$', i
            )
        ) { # e.g. 'II|III'
            chr1 <- strsplit(i, '\\|')[[1]][1]
            chr2 <- strsplit(i, '\\|')[[1]][2]
            if (!all(c(chr1, chr2) %in% seqnames(GenomeInfoDb::seqinfo(x)))) {
                stop("One or all of the provided seqnames is not found.")
            }
            si <- GenomeInfoDb::seqinfo(x)
            gr1 <- as(si[chr1], 'GRanges')
            gr2 <- as(si[chr2], 'GRanges')
            i_ <- paste(as.character(gr1), as.character(gr2), sep = '|')
            # Redirect to character (with extended i_ being like "II:xxx-xxx|III:xxx-xxx")
            return(x[i_])
        }
        else if (
            all(i %in% seqnames(GenomeInfoDb::seqinfo(x)))
        ) { # e.g. 'II'
            valid_regions_first <- re_$bin_id[as.vector(seqnames(re_)) %in% i]
            valid_regions_second <- valid_regions_first
        }
        else {
            stop("Failed to coerce i into a GRanges or GInteractions")
        }
        focus(x) <- i
    }

    else { # c('II', 'III')
        if (
            all(i %in% seqnames(GenomeInfoDb::seqinfo(x)))
        ) { 
            valid_regions_first <- re_$bin_id[as.vector(seqnames(re_)) %in% i]
            valid_regions_second <- valid_regions_first
        }
        else {
            stop("Failed to coerce i into a GRanges or GInteractions")
        }
        focus(x) <- paste(i, collapse = ', ')
    }

    sub <- ints_$bin_id1 %in% valid_regions_first & ints_$bin_id2 %in% valid_regions_second
    x[sub]
})

################################################################################
################################################################################
###############                                                  ###############
###############          METHODS FOR EXISTING GENERICS           ###############
###############                                                  ###############
################################################################################
################################################################################

#' @export
#' @rdname HiCExperiment

setMethod("fileName", "HiCExperiment", function(object) object@fileName)

#' @export
#' @rdname HiCExperiment

setMethod("interactions", signature(x = "HiCExperiment"), function(x, fillout.regions = FALSE) {
    gi <- x@interactions
    n <- names(scores(x))
    for (N in n) {
        S4Vectors::mcols(gi)[[N]] <- scores(x, N)
    }
    if (fillout.regions) {
        re <- bins(x)
        subre <- re[re$bin_id >= min(regions(gi)$bin_id) & re$bin_id <= max(regions(gi)$bin_id)]
        replaceRegions(gi) <- subre
    }
    return(gi)
})

#' @export
#' @rdname HiCExperiment

setMethod("interactions<-", signature(x = "HiCExperiment", value = "GInteractions"), function(x, value) {
    x@interactions <- value
    x
})

#' @export
#' @rdname HiCExperiment

setMethod("length", "HiCExperiment", function(x) length(interactions(x)))

#' @export
#' @rdname HiCExperiment

setReplaceMethod("$", "HiCExperiment", function(x, name, value) {
    S4Vectors::mcols(regions(interactions(x)))[, name] <- value
    return(x)
})

#' @export
#' @rdname HiCExperiment

setMethod("$", "HiCExperiment", function(x, name) {
    S4Vectors::mcols(regions(interactions(x)))[, name]
})

#' @export
#' @rdname HiCExperiment

setMethod("seqinfo", "HiCExperiment", function(x) {
    if (.is_mcool(fileName(x))) {
        si <- .cool2seqinfo(fileName(x), resolution(x))
    }
    else if (.is_cool(fileName(x))) {
        si <- .cool2seqinfo(fileName(x))
    }
    else if (.is_hic(fileName(x))) {
        si <- .hic2seqinfo(fileName(x))
    }
    else if (.is_hicpro_matrix(fileName(x)) & .is_hicpro_regions(metadata(x)$regions)) {
        si <- .hicpro2seqinfo(metadata(x)$regions)
    }
    return(si)
})

#' @export
#' @rdname HiCExperiment

setMethod("bins", "HiCExperiment", function(x) {
    if (.is_cool(fileName(x))) {
        bins <- .getCoolAnchors(
            fileName(x), resolution = NULL, balanced = TRUE
        )
        GenomeInfoDb::seqinfo(bins) <- GenomeInfoDb::seqinfo(x)
    }
    else if (.is_mcool(fileName(x))) {
        bins <- .getCoolAnchors(
            fileName(x), resolution = resolution(x), balanced = TRUE
        )
        GenomeInfoDb::seqinfo(bins) <- GenomeInfoDb::seqinfo(x)
    }
    else if (.is_hic(fileName(x))) {
        bins <- .getHicAnchors(
            fileName(x), resolution = resolution(x)
        )
        GenomeInfoDb::seqinfo(bins) <- GenomeInfoDb::seqinfo(x)
    }
    else if (.is_hicpro_matrix(fileName(x)) & .is_hicpro_regions(metadata(x)$regions)) {
        bins <- .getHicproAnchors(metadata(x)$regions)
        GenomeInfoDb::seqinfo(bins) <- GenomeInfoDb::seqinfo(x)
    }
    else {
        stop("Input format not supported")
    }
    return(bins)
})

#' @export
#' @rdname HiCExperiment

setMethod("anchors", "HiCExperiment", function(x) anchors(interactions(x)))

#' @export
#' @rdname HiCExperiment

setMethod("regions", "HiCExperiment", function(x) regions(interactions(x)))

#' @export
#' @rdname HiCExperiment

setMethod("cis", "HiCExperiment", function(x) {
    an <- anchors(x)
    sub <- GenomicRanges::seqnames(an[[1]]) == GenomicRanges::seqnames(an[[2]])
    x[which(as.vector(sub))]
})

#' @export
#' @rdname HiCExperiment

setMethod("trans", "HiCExperiment", function(x) {
    an <- anchors(x)
    sub <- GenomicRanges::seqnames(an[[1]]) != GenomicRanges::seqnames(an[[2]])
    x[which(as.vector(sub))]
})

setMethod("show", signature("HiCExperiment"), function(object) {

    if (is.null(focus(object))) {
        focus_str <- "whole genome"
    } 
    else {
        focus_str <- coords2char(focus(object))
    }
    gi <- interactions(object)
    cat(paste0(
        "`HiCExperiment` object with ", 
        format(sum(S4Vectors::mcols(gi)[,names(scores(object))[1]], na.rm = TRUE), big.mark = ","), 
        " contacts over ", 
        format(length(regions(object)), big.mark = ","), 
        " regions"
    ), '\n')
    cat('-------\n')
    cat(paste0('fileName: ', ifelse(nchar(fileName(object)) >= 1, paste0('\"', fileName(object), '\"'), "N/A")), '\n')
    cat(paste0('focus: ', ifelse(nchar(focus_str) >= 1, paste0('\"', focus_str, '\"'), "N/A")), '\n')

    ## Resolutions
    S4Vectors::coolcat("resolutions(%d): %s\n", resolutions(object))
    cat(paste0('active resolution: ', resolution(object)), '\n')

    ## Interactions
    cat(paste0('interactions: ', length(gi)), '\n')

    ## Scores
    cat(paste0('scores(', length(scores(object)), '): ', paste(names(scores(object)), collapse = " ")), '\n')

    ## topologicalFeatures
    cat(paste0('topologicalFeatures: ', paste(paste0(names(topologicalFeatures(object)), "(", lengths(topologicalFeatures(object)), ")"), collapse = " ")), '\n')

    ## Pairs
    cat(paste0('pairsFile: ', ifelse(is.null(pairsFile(object)), "N/A", pairsFile(object))), '\n')

    ## Metadata
    S4Vectors::coolcat("metadata(%d): %s\n", names(S4Vectors::metadata(object)))

})
