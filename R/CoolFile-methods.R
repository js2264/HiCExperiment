#' @title `CoolFile` methods
#' 
#' @name CoolFile-methods
#' @aliases show,CoolFile-method
#' @aliases pairsFile,CoolFile-method
#' @aliases resolution,CoolFile-method
#' 
#' @description
#' 
#' CoolFile methods.
#'
#' @param object A \code{CoolFile} object.
#' @param x A \code{CoolFile} object.
#'
#' @importFrom BiocGenerics path
#' @examples 
#' mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#' hcf <- CoolFile(mcool_path, resolution = 4000, pairs = pairs_path)
#' hcf
#' resolution(hcf)
NULL

#' @export
#' @include CoolFile-class.R 
#' @include PairsFile-class.R 

setMethod("show", signature("CoolFile"), function(object) {
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

setMethod("pairsFile", "CoolFile", function(x) {
    BiocGenerics::path(x@pairsFile)
})

#' @export

setMethod("resolution", "CoolFile", function(x) x@resolution)
