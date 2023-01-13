#' @title `AggrHiCExperiment` methods
#' 
#' @name AggrHiCExperiment
#' @rdname AggrHiCExperiment
#' @aliases slices,AggrHiCExperiment,missing-method
#' @aliases slices,AggrHiCExperiment,character-method
#' @aliases slices,AggrHiCExperiment,numeric-method
#' @aliases show,AggrHiCExperiment-method
#' 
#' @description
#' 
#' AggrHiCExperiment methods.
#' 
#' @param x,object A \code{AggrHiCExperiment} object.
#' @param name The name/index of slices to extract.
#' 
#' @include AggrHiCExperiment-class.R
#' @include HiCExperiment-methods.R
NULL

#' @export
#' @rdname AggrHiCExperiment

setMethod("slices", signature(x = "AggrHiCExperiment", name = "missing"), function(x) x@slices)

#' @export
#' @rdname AggrHiCExperiment

setMethod("slices", signature(x = "AggrHiCExperiment", name = "character"), function(x, name) {
    if (!name %in% names(slices(x))) {
        stop(paste0(name, ' not in slices.'))
    }
    return(x@slices[[name]])
})

#' @export
#' @rdname AggrHiCExperiment

setMethod("slices", signature(x = "AggrHiCExperiment", name = "numeric"), function(x, name) {
    if (name > length(slices(x))) {
        stop(paste0('Only ', length(slices(x)), ' slices in x.'))
    }
    return(x@slices[[name]])
})

#' @export
#' @rdname AggrHiCExperiment

setMethod("show", signature("AggrHiCExperiment"), function(object) {

    ntargets <- length(topologicalFeatures(object, 'targets'))
    cat(paste0('`AggrHiCExperiment` object over ', ntargets, ' targets'), '\n')
    cat('-------\n')
    cat(paste0('fileName: \"', fileName(object), '\"'), '\n')
    cat(paste0('focus: ', ntargets, ' targets'), '\n')

    ## Resolutions
    S4Vectors::coolcat("resolutions(%d): %s\n", resolutions(object))
    cat(paste0('current resolution: ', resolution(object)), '\n')

    ## Interactions
    cat(paste0('interactions: ', length(interactions(object))), '\n')

    ## Scores
    cat(paste0('scores(', length(scores(object)), '): ', paste(names(scores(object)), collapse = " ")), '\n')

    ## Slices
    cat(paste0('slices(', length(slices(object)), '): ', paste(names(slices(object)), collapse = " ")), '\n')

    ## topologicalFeatures
    cat(paste0('topologicalFeatures: ', paste(paste0(names(topologicalFeatures(object)), "(", lengths(topologicalFeatures(object)), ")"), collapse = " ")), '\n')

    ## Pairs
    cat(paste0('pairsFile: ', ifelse(is.null(pairsFile(object)), "N/A", pairsFile(object))), '\n')

    ## Metadata
    S4Vectors::coolcat("metadata(%d): %s\n", names(S4Vectors::metadata(object)))

})
