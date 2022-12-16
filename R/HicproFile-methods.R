#' @title `HicproFile` methods
#' 
#' @name HicproFile-methods
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
#' @examples 
#' hicpro_matrix_path <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_matrix')
#' hicpro_bed_path <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_bed')
#' pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#' hicpro <- HicproFile(hicpro_matrix_path, bed = hicpro_bed_path, pairs = pairs_path)
#' hicpro
#' resolution(hicpro)
#' pairsFile(hicpro)
#' S4Vectors::metadata(hicpro)
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
