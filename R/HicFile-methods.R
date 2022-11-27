#' @title `HicFile` methods
#' 
#' @name HicFile-methods
#' @aliases show,HicFile-method
#' @aliases resolution,HicFile-method
#' @aliases metadata<-,HicFile-method
#' @aliases metadata<-,HicFile,list-method
#' 
#' @description
#' 
#' HicFile methods.
#'
#' @param object A \code{HicFile} object.
#' @param x A \code{HicFile} object.
#'
#' @importFrom BiocGenerics path
#' @examples 
#' hic_path <- HiContactsData::HiContactsData('yeast_wt', 'hic')
#' pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#' hcf <- HicFile(hic_path, resolution = 4000, pairs = pairs_path)
#' hcf
#' resolution(hcf)
NULL

#' @export
#' @include HicFile-class.R 

setMethod("show", signature("HicFile"), function(object) {
    r <- BiocIO::resource(object)
    res <- resolution(object)
    if (is.null(res)) res = lsHicResolutions(r)[1]
    if (!S4Vectors::isSingleString(r))
        stop('"filename" must be a single string, specifiying a path')
    cat(class(object), "object\n.hic file:", r, '\n') 
    cat("resolution:", res, '\n')
    S4Vectors::coolcat("metadata(%d): %s\n", names(S4Vectors::metadata(object)))
})

#' @export

setMethod("resolution", "HicFile", function(x) x@resolution)

#' @export

setMethod("metadata<-", signature(x = "HicFile", value = "list"), function(x, value) {
    x@metadata <- value
    x
})
