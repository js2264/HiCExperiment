#' @title `HiCoolFile` methods
#' 
#' @name HiCoolFile-methods
#' @aliases show,HiCoolFile-method
#' @aliases pairsFile,HiCoolFile-method
#' 
#' @description
#' 
#' HiCoolFile methods.
#'
#' @param object A \code{HiCoolFile} object.
#' @param x A \code{HiCoolFile} object.
#'
#' @importFrom BiocGenerics path
#' @examples 
#' mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#' hcf <- HiCoolFile(mcool_path, resolution = 4000, pairs = pairs_path)
#' hcf
#' resolution(hcf)
NULL

#' @export

setMethod("show", signature("HiCoolFile"), function(object) {
    r <- BiocIO::resource(object)
    res <- resolution(object)
    if (is.null(res)) res = lsCoolResolutions(r)[1]
    if (!S4Vectors::isSingleString(r))
        stop('"filename" must be a single string, specifiying a path')
    cat(
        class(object), "object\nmcool file:", r, 
        "\nresolution:", res, 
        "\npairs file:", pairsFile(object), 
        "\nlog:", S4Vectors::metadata(object)$log, 
        "\nstats:", S4Vectors::metadata(object)$stats, '\n'  
    )
})

#' @export

setMethod("pairsFile", "HiCoolFile", function(x) {
    BiocGenerics::path(x@pairsFile)
})
