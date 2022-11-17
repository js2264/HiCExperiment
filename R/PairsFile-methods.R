#' @title `pairsFile` methods
#' 
#' @name pairsFile-methods
#' @aliases pairsFile,PairsFile-method
#' 
#' @description
#' 
#' PairsFile methods
#' 
#' @param x Path to a pairs file
#'
#' @importFrom BiocGenerics path
#' @examples 
#' pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#' pf <- PairsFile(pairs_path)
#' pairsFile(pf)
NULL

#' @export

setMethod("pairsFile", "PairsFile", function(x) {
    BiocGenerics::path(x)
})
