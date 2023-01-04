#' @title `HicproFile` methods
#' 
#' @name HicproFile-methods
#' @rdname HicproFile-class
#' @aliases show,HicproFile-method
#' 
#' @description
#' 
#' HicproFile methods.
#'
#' @param object A \code{HicproFile} object.
#'
#' @importFrom BiocGenerics path
#' @include HicproFile-class.R 
#' @include PairsFile-class.R 
NULL

#' @export

setMethod("show", signature("HicproFile"), function(object) {
    r <- BiocIO::resource(object)
    res <- resolution(object)
    cat(class(object), "object\nHiC-Pro files:\n  $ matrix:  ", r, '\n  $ regions: ', object@bed, '\n') 
    cat("resolution:", res, '\n')
    cat("pairs file:", pairsFile(object), '\n')
    S4Vectors::coolcat("metadata(%d): %s\n", names(S4Vectors::metadata(object)))
})
