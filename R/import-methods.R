#' @title HiCExperiment import methods
#' 
#' @name import-methods
#' @aliases import
#' @aliases import,CoolFile-method
#' @aliases import,PairsFile-method
#' @aliases import,CoolFile,ANY,ANY-method
#' @aliases import,PairsFile,ANY,ANY-method
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
#' @param format The format of the output. If missing and 'con' is a filename, 
#'    the format is derived from the file extension. 
#'    This argument is unnecessary when 'con' is a derivative of 'BiocFile'.
#' @param text If 'con' is missing, this can be a character vector directly 
#'    providing the string data to import.
#' @param ... e.g. `resolution = ...`; parameters to pass to the 
#'    format-specific method.
#' 
#' @usage import(con, format, text, ...)
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
