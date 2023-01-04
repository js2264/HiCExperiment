#' @title `pairsFile` methods
#' 
#' @name PairsFile-class
#' @rdname PairsFile-class
#' @aliases pairsFile,PairsFile-method
#' 
#' @description
#' 
#' PairsFile methods
#' 
#' @param x Path to a pairs file
#'
#' @importFrom BiocGenerics path
NULL

#' @export

setMethod("pairsFile", "PairsFile", function(x) {
    BiocGenerics::path(x)
})
