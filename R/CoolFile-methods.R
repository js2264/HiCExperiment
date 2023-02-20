#' @title `CoolFile` methods
#' 
#' @name CoolFile-methods
#' @rdname CoolFile-class
#' @aliases show,CoolFile-method
#' 
#' @description
#' 
#' CoolFile methods.
#'
#' @param object A \code{CoolFile} object.
#'
#' @importFrom S4Vectors metadata
#' @importFrom BiocGenerics path
#' @include CoolFile-class.R 
#' @include PairsFile-class.R 
NULL

#' @export

setMethod("show", signature("CoolFile"), function(object) {
    r <- BiocIO::resource(object)
    res <- resolution(object)
    if (is.null(res)) res = .lsCoolResolutions(r)[1]
    if (!S4Vectors::isSingleString(r))
        stop('"filename" must be a single string, specifiying a path')
    cat(class(object), "object\n.mcool file:", r, '\n') 
    cat("resolution:", res, '\n')
    cat("pairs file:", pairsFile(object), '\n')
    S4Vectors::coolcat("metadata(%d): %s\n", names(S4Vectors::metadata(object)))
})
