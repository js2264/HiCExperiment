#' @title `CoolFile` methods
#' 
#' @name CoolFile-methods
#' @aliases show,CoolFile-method
#' @aliases resolution,CoolFile-method
#' 
#' @description
#' 
#' CoolFile methods.
#'
#' @param object A \code{CoolFile} object.
#' @param x A \code{CoolFile} object.
#'
#' @examples 
#' cool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' cf <- CoolFile(cool_path, resolution = 4000)
#' cf
#' resolution(cf)
NULL

#' @export

setMethod("show", signature("CoolFile"), function(object) {
    r <- BiocIO::resource(object)
    res <- resolution(object)
    if (is.null(res)) res = lsCoolResolutions(r)[1]
    if (!S4Vectors::isSingleString(r))
        stop('"filename" must be a single string, specifiying a path')
    cat(class(object), "object\nresource:", r, "\nresolution:", res, "\n")
})

#' @export
setMethod("resolution", "CoolFile", function(x) x@resolution)
