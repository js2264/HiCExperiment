#' @title `AggrHiCExperiment` S4 class
#' 
#' @name AggrHiCExperiment
#' @rdname AggrHiCExperiment
#' 
#' @description
#' 
#' The `AggrHiCExperiment` extends `HiCExperiment` class
#'
#' @slot fileName Path of Hi-C contact file
#' @slot resolutions Resolutions available in the Hi-C contact file.
#' @slot resolution Current resolution
#' @slot interactions Genomic Interactions extracted from the Hi-C contact file
#' @slot scores Available interaction scores. 
#' @slot topologicalFeatures Topological features associated with the dataset 
#'   (e.g. loops (\<Pairs\>), borders (\<GRanges\>), 
#'   viewpoints (\<GRanges\>), etc...)
#' @slot pairsFile Path to the .pairs file associated with the Hi-C contact file
#' @slot metadata metadata associated with the Hi-C contact file.
#' 
#' @include AggrHiCExperiment-class.R
#' @include HiCExperiment-methods.R
#' 
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

    nsnippets <- length(topologicalFeatures(object, 'snippets'))
    cat(glue::glue('`AggrHiCExperiment` object over {nsnippets} snippets'), '\n')
    cat('-------\n')
    cat(glue::glue('fileName: "{fileName(object)}"'), '\n')
    cat(glue::glue('focus: {nsnippets} snippets'), '\n')

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
