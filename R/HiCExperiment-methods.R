#' @title `HiCExperiment` methods
#' 
#' @name HiCExperiment-methods
#' @aliases resolutions,HiCExperiment-method
#' @aliases resolution,HiCExperiment-method
#' @aliases focus,HiCExperiment-method
#' @aliases focus<-,HiCExperiment-method
#' @aliases focus<-,HiCExperiment,character-method
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
#' @aliases [,HiCExperiment,logical,ANY,ANY-method
#' @aliases [,HiCExperiment,character,ANY,ANY-method
#' @aliases seqinfo,HiCExperiment-method
#' @aliases bins,HiCExperiment-method
#' @aliases anchors,HiCExperiment-method
#' @aliases regions,HiCExperiment-method
#' @aliases show,HiCExperiment-method
#' @aliases coerce,HiCExperiment,GInteractions-method
#' @aliases coerce,HiCExperiment,ContactMatrix-method
#' @aliases coerce,HiCExperiment,matrix-method
#' @aliases coerce,HiCExperiment,data.frame-method
#' @aliases as.matrix,HiCExperiment-method
#' @aliases as.data.frame,HiCExperiment-method
#' 
#' @description
#' 
#' HiCExperiment methods.
#'
#' @param x A \code{HiCExperiment} object.
#' @param object A \code{HiCExperiment} object.
#' @param name ...
#' @param value ...
#' @param i a range or boolean vector.
#' 
#' @importMethodsFrom BiocGenerics fileName
#' @importFrom InteractionSet reduceRegions
#' @importFrom GenomeInfoDb seqinfo
#' 
#' @examples 
#' library(HiCExperiment)
#' mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' contacts_yeast <- HiCExperiment(mcool_path, resolution = 16000)
#' resolutions(contacts_yeast)
#' resolution(contacts_yeast)
#' focus(contacts_yeast)
#' scores(contacts_yeast)
#' tail(scores(contacts_yeast, 1))
#' tail(scores(contacts_yeast, 'balanced'))
#' scores(contacts_yeast, 'test') <- runif(length(contacts_yeast))
#' tail(scores(contacts_yeast, 'test'))
#' full_contacts_yeast <- full_contacts_yeast()
#' topologicalFeatures(full_contacts_yeast)
#' topologicalFeatures(full_contacts_yeast, 1)
#' topologicalFeatures(full_contacts_yeast, 'centromeres')
#' data(centros_yeast)
#' topologicalFeatures(contacts_yeast, 'centromeres') <- centros_yeast
#' topologicalFeatures(contacts_yeast, 'centromeres')
#' pairsFile(full_contacts_yeast)
#' fileName(contacts_yeast)
#' interactions(contacts_yeast)
#' length(contacts_yeast)
#' contacts_yeast[seq_len(10)]
#' seqinfo(contacts_yeast)
#' bins(contacts_yeast)
#' anchors(contacts_yeast)
#' regions(contacts_yeast)
#' contacts_yeast
#' as(contacts_yeast, 'GInteractions')
#' as(contacts_yeast, 'ContactMatrix')
#' as(contacts_yeast, 'matrix')[seq_len(10), seq_len(10)]
#' as(contacts_yeast, 'data.frame')[seq_len(10), seq_len(10)]
NULL

################################################################################
################################################################################
###############                                                  ###############
###############          METHODS FOR NEW GENERICS                ###############
###############                                                  ###############
################################################################################
################################################################################

#' @export

setMethod("resolutions", "HiCExperiment", function(x) x@resolutions)

#' @export

setMethod("resolution", "HiCExperiment", function(x) x@resolution)

#' @export

setMethod("focus", "HiCExperiment", function(x) x@focus)

setMethod("focus<-", signature(x = "HiCExperiment", value = "character"), function(x, value) {
    x@focus <- value
    x
})

#' @export

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

#' @export

setMethod("scores<-", c(x = "HiCExperiment", name = "character", value = "numeric"), function(x, name, value) {
    x@scores[[name]] <- value
    return(x)
})

#' @export

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

#' @export

setMethod("topologicalFeatures<-", signature(x = "HiCExperiment", name = "character", value = "GRangesOrGInteractions"), function(x, name, value) {
    x@topologicalFeatures[[name]] <- value
    return(x)
})

#' @export

setMethod("pairsFile", "HiCExperiment", function(x) {
    x@pairsFile
})

#' @export

setMethod("pairsFile<-", signature(x = "HiCExperiment", value = "character"), function(x, value) {
    if (!file.exists(value)) {
        stop("Provided pairsFile does not exist. Aborting now.")
    }
    x@pairsFile <- value
    x
})

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

#' @export

setMethod("fileName", "HiCExperiment", function(object) object@fileName)

#' @export

setMethod("interactions", "HiCExperiment", function(x) x@interactions)

setMethod("interactions<-", signature(x = "HiCExperiment", value = "GInteractions"), function(x, value) {
    x@interactions <- value
    x
})

#' @export

setMethod("length", "HiCExperiment", function(x) length(interactions(x)))

#' @export

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

#' @export

setMethod("seqinfo", "HiCExperiment", function(x) {
    if (is_mcool(fileName(x))) {
        si <- cool2seqinfo(fileName(x), resolution(x))
    }
    else {
        si <- cool2seqinfo(fileName(x))
    }
    return(si)
})

#' @export

setMethod("bins", "HiCExperiment", function(x) {
    bins <- getAnchors(
        fileName(x), resolution = resolution(x), balanced = FALSE
    )
    GenomeInfoDb::seqinfo(bins) <- GenomeInfoDb::seqinfo(x)
    return(bins)
})

#' @export

setMethod("anchors", "HiCExperiment", function(x) anchors(interactions(x)))

#' @export

setMethod("regions", "HiCExperiment", function(x) regions(interactions(x)))

#' @export

setMethod("show", signature("HiCExperiment"), function(object) {

    if (is.null(focus(object))) {
        focus_str <- "whole genome"
    } 
    else {
        focus_str <- coords2char(focus(object))
    }

    cat(glue::glue(
        '`HiCExperiment` object with {format(length(interactions(object)), big.mark = ",")} interactions over {format(length(regions(object)), big.mark = ",")} regions'
    ), '\n')
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

#' @export

setAs("HiCExperiment", "GInteractions", function(from) interactions(from))

#' @export

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

#' @export

setAs("HiCExperiment", "matrix", function(from) {
    as(from, "ContactMatrix") |> cm2matrix()
})

#' @export

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

setMethod("as.matrix", "HiCExperiment", function(x) {
    as(x, 'matrix')
})

#' @export

setMethod("as.data.frame", "HiCExperiment", function(x) {
    as(x, 'data.frame')
})