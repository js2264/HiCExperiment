#' @title HiCExperiment import methods
#' 
#' @name import-methods
#' @aliases import
#' @aliases import,CoolFile-method
#' @aliases import,HiCoolFile-method
#' @aliases import,PairsFile-method
#' 
#' @description 
#' 
#' Import methods for data structures implemented in the HiCExperiment package. 
#' HiCExperiment package implements methods to faciliate 
#' the import of (m)cool files and pairs files in R, as HiCExperiment  
#' or GenomicInteractions objects.
#' 
#' @param con Path or connection to a cool, mcool or pairs file. Con argument 
#'   can also be a CoolFile or PairsFile object. 
#' @param ... e.g. `resolution = ...`; parameters to pass to the 
#'    format-specific method.
#' 
#' @usage import(con, ...)
#' 
#' @importFrom BiocIO import
#' @importFrom BiocGenerics path
#' @exportMethod import 
NULL

#' @export

setMethod('import', 'CoolFile', function(con, ...) {

    path <- BiocGenerics::path(con)
    stopifnot(file.exists(path))
    params <- list(...)
    if ('resolution' %in% names(params)) {
        check_cool_format(path, params[['resolution']])
        HiCExperiment(path, ...)
    }
    else {
        check_cool_format(path, resolution(con))
        HiCExperiment(con, ...)
    }

})

#' @export

setMethod('import', 'HiCoolFile', function(con, ...) {

    path <- BiocGenerics::path(con)
    stopifnot(file.exists(path))
    params <- list(...)
    if ('resolution' %in% names(params)) {
        check_cool_format(path, params[['resolution']])
        HiCExperiment(
            path, ..., 
            pairsFile = pairsFile(con), 
            metadata = S4Vectors::metadata(con)
        )
    }
    else {
        check_cool_format(path, resolution(con))
        HiCExperiment(
            con, ..., 
            pairsFile = pairsFile(con), 
            metadata = S4Vectors::metadata(con)
        )
    }

})

#' @export

setMethod('import', signature = 'PairsFile', function(con, ...) {

    con <- BiocGenerics::path(con)
    stopifnot(file.exists(con))
    pairs2gi(con, ...)

})
