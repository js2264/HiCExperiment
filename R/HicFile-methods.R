#' @title `HicFile` methods
#' 
#' @name HicFile-methods
#' @aliases show,HicFile-method
#' 
#' @description
#' 
#' HicFile methods.
#'
#' @param object A \code{HicFile} object.
#'
#' @importFrom BiocGenerics path
#' @include HicFile-class.R 
#' @include PairsFile-class.R 
#' @examples 
#' hic_path <- HiContactsData::HiContactsData('yeast_wt', 'hic')
#' pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#' hic <- HicFile(
#'   hic_path, resolution = 4000, 
#'   pairsFile = pairs_path, metadata = list(info = 'Processed with Juicer')
#' )
#' hic
#' resolution(hic)
#' pairsFile(hic)
#' metadata(hic)
NULL

#' @export

setMethod("show", signature("HicFile"), function(object) {
    r <- BiocIO::resource(object)
    res <- resolution(object)
    if (is.null(res)) res = lsHicResolutions(r)[1]
    if (!S4Vectors::isSingleString(r))
        stop('"filename" must be a single string, specifiying a path')
    cat(class(object), "object\n.hic file:", r, '\n') 
    cat("resolution:", res, '\n')
    cat("pairs file:", pairsFile(object), '\n')
    S4Vectors::coolcat("metadata(%d): %s\n", names(S4Vectors::metadata(object)))
})
