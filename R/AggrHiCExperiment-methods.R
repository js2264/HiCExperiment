#' @title `AggrHiCExperiment` methods
#' 
#' @name AggrHiCExperiment-methods
#' @aliases slices,AggrHiCExperiment,missing-method
#' @aliases slices,AggrHiCExperiment,character-method
#' @aliases slices,AggrHiCExperiment,numeric-method
#' @aliases show,AggrHiCExperiment-method
#' 
#' @description
#' 
#' AggrHiCExperiment methods.
#' 
#' @param x A \code{AggrHiCExperiment} object.
#' @param name The name/index of slices to extract.
#' 
#' @include AggrHiCExperiment-class.R
#' @include HiCExperiment-methods.R
NULL

#' @export

setMethod("slices", signature(x = "AggrHiCExperiment", name = "missing"), function(x) x@slices)
setMethod("slices", signature(x = "AggrHiCExperiment", name = "character"), function(x, name) {
    if (!name %in% names(slices(x))) {
        stop(paste0(name, ' not in slices.'))
    }
    return(x@slices[[name]])
})
setMethod("slices", signature(x = "AggrHiCExperiment", name = "numeric"), function(x, name) {
    if (name > length(slices(x))) {
        stop(paste0('Only ', length(slices(x)), ' slices in x.'))
    }
    return(x@slices[[name]])
})

#' @export

setMethod("show", signature("AggrHiCExperiment"), function(object) {

    ntargets <- length(topologicalFeatures(object, 'targets'))
    cat(glue::glue('`AggrHiCExperiment` object over {ntargets} targets'), '\n')
    cat('-------\n')
    cat(glue::glue('fileName: "{fileName(object)}"'), '\n')
    cat(glue::glue('focus: {ntargets} targets'), '\n')

    ## Resolutions
    S4Vectors::coolcat("resolutions(%d): %s\n", resolutions(object))
    cat(glue::glue('current resolution: {resolution(object)}'), '\n')

    ## Interactions
    cat(glue::glue('interactions: {length(interactions(object))}'), '\n')

    ## Scores
    cat(glue::glue('scores({length(scores(object))}): {paste(names(scores(object)), collapse = " ")}'), '\n')

    ## Slices
    cat(glue::glue('slices({length(slices(object))}): {paste(names(slices(object)), collapse = " ")}'), '\n')

    ## topologicalFeatures
    cat(glue::glue('topologicalFeatures: {paste(paste0(names(topologicalFeatures(object)), "(", lengths(topologicalFeatures(object)), ")"), collapse = " ")}'), '\n')

    ## Pairs
    cat(glue::glue('pairsFile: {ifelse(is.null(pairsFile(object)), "N/A", pairsFile(object))}'), '\n')

    ## Metadata
    S4Vectors::coolcat("metadata(%d): %s\n", names(S4Vectors::metadata(object)))

})
