#' @rdname PairsFile
#'
#' @name pairsFile
#' @docType methods
#' @aliases pairsFile,PairsFile-method
#'
#' @param x 
#'
#' @importFrom BiocGenerics path
#' @export

setMethod("pairsFile", signature("PairsFile", "missing"), function(x) {
    BiocGenerics::path(x)
})

#' @importFrom BiocGenerics path
#' @export

setMethod("pairsFile", signature("HiCoolFile", "missing"), function(x) {
    BiocGenerics::path(x@pairsFile)
})
