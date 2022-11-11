################################################################################
################################################################################
###############                                                  ###############
###############          METHODS FOR NEW GENERICS                ###############
###############                                                  ###############
################################################################################
################################################################################

#' @rdname HiCExperiment
#'
#' @name resolutions
#' @docType methods
#' @aliases resolutions,HiCExperiment-method
#'
#' @param x A \code{HiCExperiment} object.
#' 
#' @export
#' @examples 
#' library(HiCExperiment)
#' contacts_yeast <- contacts_yeast()
#' resolutions(contacts_yeast)

setMethod("resolutions", "HiCExperiment", function(x) x@resolutions)

#' @rdname HiCExperiment
#'
#' @name resolution
#' @docType methods
#' @aliases resolution,HiCExperiment-method
#'
#' @export
#' @examples 
#' resolution(contacts_yeast)

setMethod("resolution", "HiCExperiment", function(x) x@resolution)

#' @rdname HiCExperiment
#'
#' @name focus
#' @docType methods
#' @aliases focus,HiCExperiment-method
#'
#' @export
#' @examples 
#' focus(contacts_yeast)

setMethod("focus", "HiCExperiment", function(x) x@focus)

#' @rdname HiCExperiment
#'
#' @name focus<-
#' @docType methods
#' @aliases focus<-,HiCExperiment,character-method
#'
#' @param name name
#' @param value value
#'
#' @export

setMethod("focus<-", signature(x = "HiCExperiment", value = "character"), function(x, value) {
    x@focus <- value
    x
})

#' @rdname HiCExperiment
#'
#' @name scores
#' @docType methods
#' @aliases scores,HiCExperiment,missing-method
#' @aliases scores,HiCExperiment,character-method
#' @aliases scores,HiCExperiment,numeric-method
#'
#' @export
#' @examples 
#' scores(contacts_yeast)
#' tail(scores(contacts_yeast, 1))
#' tail(scores(contacts_yeast, 'balanced'))

setMethod("scores", signature(x = "HiCExperiment", name = "missing"), function(x) x@scores)
setMethod("scores", signature(x = "HiCExperiment", name = "character"), function(x, name) {
    if (!name %in% names(scores(x))) {
        stop(paste0(name, ' not in scores.'))
    }
    return(x@scores[[name]])
})
setMethod("scores", signature(x = "HiCExperiment", name = "numeric"), function(x, name) {
    if (name > length(scores(x))) {
        stop(paste0('Only ', length(scores(x)), ' scores in x.'))
    }
    return(x@scores[[name]])
})

#' @rdname HiCExperiment
#'
#' @name scores<-
#' @docType methods
#' @aliases scores<-,HiCExperiment,character,numeric-method
#'
#' @param name name
#' @param value value
#'
#' @export
#' @examples 
#' scores(contacts_yeast, 'test') <- runif(length(contacts_yeast))
#' tail(scores(contacts_yeast, 'test'))

setMethod("scores<-", c(x = "HiCExperiment", name = "character", value = "numeric"), function(x, name, value) {
    x@scores[[name]] <- value
    return(x)
})

#' @rdname HiCExperiment
#' 
#' @name topologicalFeatures
#' @docType methods
#' @aliases topologicalFeatures,HiCExperiment,missing-method
#' @aliases topologicalFeatures,HiCExperiment,character-method
#' @aliases topologicalFeatures,HiCExperiment,numeric-method
#'
#' @export
#' @examples 
#' full_contacts_yeast <- full_contacts_yeast()
#' topologicalFeatures(full_contacts_yeast)
#' topologicalFeatures(full_contacts_yeast, 1)
#' topologicalFeatures(full_contacts_yeast, 'centromeres')

setMethod("topologicalFeatures", signature(x = "HiCExperiment", name = "missing"), function(x) {
    S4Vectors::SimpleList(as.list(x@topologicalFeatures))
})
setMethod("topologicalFeatures", signature(x = "HiCExperiment", name = "character"), function(x, name) {
    if (!name %in% names(topologicalFeatures(x))) {
        stop(paste0(name, ' not in topologicalFeatures.'))
    }
    x@topologicalFeatures[[name]]
})
setMethod("topologicalFeatures", signature(x = "HiCExperiment", name = "numeric"), function(x, name) {
    if (name > length(topologicalFeatures(x))) {
        stop(paste0('Only ', length(topologicalFeatures(x)), ' topologicalFeatures in x.'))
    }
    x@topologicalFeatures[[name]]
})

#' @rdname HiCExperiment
#' 
#' @name topologicalFeatures<-
#' @docType methods
#' @aliases topologicalFeatures<-,HiCExperiment,character,GRangesOrGInteractions-method
#'
#' @param name name
#' @param value value
#'
#' @export
#' @examples 
#' data(centros_yeast)
#' topologicalFeatures(contacts_yeast, 'centromeres') <- centros_yeast
#' topologicalFeatures(contacts_yeast, 'centromeres')

setMethod("topologicalFeatures<-", signature(x = "HiCExperiment", name = "character", value = "GRangesOrGInteractions"), function(x, name, value) {
    x@topologicalFeatures[[name]] <- value
    return(x)
})

#' @rdname HiCExperiment
#' 
#' @name pairsFile
#' @docType methods
#' @aliases pairsFile,HiCExperiment-method
#'
#' @export
#' @examples 
#' pairsFile(full_contacts_yeast)

setMethod("pairsFile", "HiCExperiment", function(x) {
    x@pairsFile
})

#' @rdname HiCExperiment
#' 
#' @name pairsFile<-
#' @docType methods
#' @aliases pairsFile<-,HiCExperiment,character-method
#'
#' @param name name
#' @param value value
#'
#' @export

setMethod("pairsFile<-", signature(x = "HiCExperiment", value = "character"), function(x, value) {
    if (!file.exists(value)) {
        stop("Provided pairsFile does not exist. Aborting now.")
    }
    x@pairsFile <- value
    x
})

#' @rdname HiCExperiment
#' 
#' @name metadata<-
#' @docType methods
#' @aliases metadata<-,HiCExperiment,list-method
#'
#' @param name name
#' @param value value
#'
#' @export

setMethod("metadata<-", signature(x = "HiCExperiment", value = "list"), function(x, value) {
    x@metadata <- value
    x
})

################################################################################
################################################################################
###############                                                  ###############
###############          METHODS FOR EXISTING GENERICS           ###############
###############                                                  ###############
################################################################################
################################################################################

#' @rdname HiCExperiment
#'
#' @name fileName
#' @docType methods
#' @aliases fileName,HiCExperiment-method
#'
#' @importMethodsFrom BiocGenerics fileName
#' @export
#' @examples 
#' fileName(contacts_yeast)

setMethod("fileName", "HiCExperiment", function(object) object@fileName)

#' @rdname HiCExperiment
#'
#' @name interactions
#' @docType methods
#' @aliases interactions,HiCExperiment-method
#'
#' @export
#' @examples 
#' interactions(contacts_yeast)

setMethod("interactions", "HiCExperiment", function(x) x@interactions)

#' @rdname HiCExperiment
#'
#' @name interactions<-
#' @docType methods
#' @aliases interactions<-,HiCExperiment,GInteractions-method
#'
#' @param name name
#' @param value value
#'
#' @export

setMethod("interactions<-", signature(x = "HiCExperiment", value = "GInteractions"), function(x, value) {
    x@interactions <- value
    x
})

#' @rdname HiCExperiment
#'
#' @name length
#' @docType methods
#' @aliases length,HiCExperiment-method
#'
#' @export
#' @examples 
#' length(contacts_yeast)

setMethod("length", "HiCExperiment", function(x) length(interactions(x)))

#' @rdname HiCExperiment
#'
#' @name [
#' @docType methods
#' @aliases [,HiCExperiment,numeric,ANY,ANY-method
#' @aliases [,HiCExperiment,logical,ANY,ANY-method
#' @aliases [,HiCExperiment,character,ANY,ANY-method
#'
#' @param i a range or boolean vector.
#'
#' @importFrom InteractionSet reduceRegions
#' @importFrom GenomeInfoDb seqinfo
#' @export
#' @examples 
#' contacts_yeast[seq_len(10)]

setMethod("[", signature("HiCExperiment", "numeric"), function(x, i) {
    interactions(x) <- InteractionSet::reduceRegions(
        interactions(x)[i]
    )
    for (n in names(scores(x))) {
        scores(x, n) <- scores(x, n)[i]
    }
    return(x)
})
setMethod("[", signature("HiCExperiment", "logical"), function(x, i) {
    interactions(x) <- InteractionSet::reduceRegions(
        interactions(x)[i]
    )
    for (n in names(scores(x))) {
        scores(x, n) <- scores(x, n)[i]
    }
    return(x)
})
setMethod("[", signature("HiCExperiment", "character"), function(x, i) {
    re_ <- regions(x)
    ints_ <- interactions(x)
    if (length(i) == 1) {
        if (grepl(
            '[A-Za-z0-9]*:[0-9]*-[0-9]* [xX/-;] [A-Za-z0-9]*:[0-9]*-[0-9]*$', i
        )) {
            i_ <- char2coords(i)
            valid_regions_first <- subsetByOverlaps(
                re_, S4Vectors::first(i_), type = 'within'
            )$bin_id
            valid_regions_second <- subsetByOverlaps(
                re_, S4Vectors::second(i_), type = 'within'
            )$bin_id
        }
        else if (grepl(
            '[A-Za-z0-9]*:[0-9]*-[0-9]*$', i
        )) {
            i_ <- char2coords(i)
            valid_regions_first <- subsetByOverlaps(
                re_, S4Vectors::first(i_), type = 'within'
            )$bin_id
            valid_regions_second <- valid_regions_first
        }
        else if (
            i %in% seqnames(GenomeInfoDb::seqinfo(x))
        ){
            valid_regions_first <- re_$bin_id[as.vector(seqnames(re_)) %in% i]
            valid_regions_second <- valid_regions_first
        }
        else {
            stop("Failed to coerce i into a Pairs/GRanges/chr.")
        }
    }
    else {
        if (
            all(i %in% seqnames(GenomeInfoDb::seqinfo(x)))
        ){
            valid_regions_first <- re_$bin_id[as.vector(seqnames(re_)) %in% i]
            valid_regions_second <- valid_regions_first
        }
        else {
            stop("Failed to coerce i into a valid Pairs/GRanges/chr.")
        }
    }
    sub <- ints_$bin_id1 %in% valid_regions_first & ints_$bin_id2 %in% valid_regions_second
    interactions(x) <- InteractionSet::reduceRegions(
        ints_[sub]
    )
    for (n in names(scores(x))) {
        scores(x, n) <- scores(x, n)[sub]
    }
    focus(x) <- i
    return(x)
})

#' @rdname HiCExperiment
#'
#' @name seqinfo,HiCExperiment-method
#' @docType methods
#'
#' @export
#' @examples 
#' seqinfo(contacts_yeast)

setMethod("seqinfo", "HiCExperiment", function(x) {
    if (is_mcool(fileName(x))) {
        si <- cool2seqinfo(fileName(x), resolution(x))
    }
    else {
        si <- cool2seqinfo(fileName(x))
    }
    return(si)
})

#' @rdname HiCExperiment
#'
#' @name bins
#' @docType methods
#' @aliases bins,HiCExperiment-method
#'
#' @importFrom GenomeInfoDb seqinfo
#' @export
#' @examples 
#' bins(contacts_yeast)

setMethod("bins", "HiCExperiment", function(x) {
    bins <- getAnchors(
        fileName(x), resolution = resolution(x), balanced = FALSE
    )
    GenomeInfoDb::seqinfo(bins) <- GenomeInfoDb::seqinfo(x)
    return(bins)
})

#' @rdname HiCExperiment
#' 
#' @name anchors
#' @docType methods
#' @aliases anchors,HiCExperiment-method
#'
#' @export
#' @examples 
#' anchors(contacts_yeast)

setMethod("anchors", "HiCExperiment", function(x) anchors(interactions(x)))

#' @rdname HiCExperiment
#' 
#' @name regions
#' @docType methods
#' @aliases regions,HiCExperiment-method
#'
#' @export
#' @examples 
#' regions(contacts_yeast)

setMethod("regions", "HiCExperiment", function(x) regions(interactions(x)))

#' @rdname HiCExperiment
#' 
#' @name summary
#' @docType methods
#' @aliases summary,HiCExperiment-method
#'
#' @export
#' @examples 
#' summary(contacts_yeast)

setMethod("summary", "HiCExperiment", function(object) {
    cat(glue::glue(
        '`HiCExperiment` object with {format(length(interactions(object)), big.mark = ",")} interactions over {format(length(regions(object)), big.mark = ",")} regions'
    ), '\n')
})

#' @rdname HiCExperiment
#' 
#' @name show
#' @docType methods
#' @aliases show,HiCExperiment-method
#'
#' @param object A \code{HiCExperiment} object.
#'
#' @export
#' @examples 
#' show(contacts_yeast)

setMethod("show", signature("HiCExperiment"), function(object) {

    if (is.null(focus(object))) {
        focus_str <- "whole genome"
    } 
    else {
        focus_str <- coords2char(focus(object))
    }

    cat(summary(object))
    cat('-------\n')
    cat(glue::glue('fileName: "{fileName(object)}"'), '\n')
    cat(glue::glue('focus: "{focus_str}"'), '\n')

    ## Resolutions
    S4Vectors::coolcat("resolutions(%d): %s\n", resolutions(object))
    cat(glue::glue('current resolution: {resolution(object)}'), '\n')

    ## Interactions
    cat(glue::glue('interactions: {length(interactions(object))}'), '\n')

    ## Scores
    cat(glue::glue('scores({length(scores(object))}): {paste(names(scores(object)), collapse = " ")}'), '\n')

    ## topologicalFeatures
    cat(glue::glue('topologicalFeatures: {paste(paste0(names(topologicalFeatures(object)), "(", lengths(topologicalFeatures(object)), ")"), collapse = " ")}'), '\n')

    ## Pairs
    cat(glue::glue('pairsFile: {ifelse(is.null(pairsFile(object)), "N/A", pairsFile(object))}'), '\n')

    ## Metadata
    S4Vectors::coolcat("metadata(%d): %s\n", names(S4Vectors::metadata(object)))

})

################################################################################
################################################################################
###############                                                  ###############
###############                    COERCING                      ###############
###############                                                  ###############
################################################################################
################################################################################

#' @rdname HiCExperiment
#' 
#' @name setAs
#' @docType methods
#' @aliases setAs,HiCExperiment-method
#'
#' @export
#' @examples 
#' as(contacts_yeast, 'GInteractions')
#' as(contacts_yeast, 'ContactMatrix')
#' as(contacts_yeast, 'matrix')[seq_len(10), seq_len(10)]
#' as(contacts_yeast, 'data.frame')[seq_len(10), seq_len(10)]

setAs("HiCExperiment", "GInteractions", function(from) interactions(from))
setAs("HiCExperiment", "ContactMatrix", function(from) {
    if ('balanced' %in% names(scores(from))) {
        x <- interactions(from)
        x$score <- scores(from, 'balanced')
        gi2cm(x)
    } 
    else {
        x <- interactions(from)
        x$score <- scores(from, 1)
        gi2cm(x)
    }
})
setAs("HiCExperiment", "matrix", function(from) {
    as(from, "ContactMatrix") |> cm2matrix()
})
setAs("HiCExperiment", "data.frame", function(from) {
    x <- interactions(from)
    x <- as.data.frame(x)
    x <- x[, !colnames(x) %in% c("chr1", "chr2", "bin_id1.1", "bin_id2.1")]
    for (n in names(scores(from))) {
        x[[n]] <- scores(from, n)
    }
    return(x)
})
