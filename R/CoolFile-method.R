#' @rdname CoolFile
#'
#' @name resolution
#' @docType methods
#' @aliases resolution,CoolFile-method
#'
#' @export
#' @examples 
#' cool_path <- HiContactsData::HiContactsData('yeast_wt', 'cool')
#' cfile <- CoolFile(cool_path)
#' resolution(cfile)
#' mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' cfile <- CoolFile(mcool_path, resolution)
#' resolution(cfile)

setMethod("resolution", "CoolFile", function(x) x@resolution)

#' @rdname CoolFile
#'
#' @name show
#' @docType methods
#' @aliases show,CoolFile-method
#'
#' @param object A \code{CoolFile} object.
#'
#' @export
#' @examples 
#' show(contacts_yeast)

setMethod("show", signature("CoolFile"), function(object) {
    r <- BiocIO::resource(object)
    res <- resolution(object)
    if (is.null(res)) res = lsCoolResolutions(r)[1]
    if (!S4Vectors::isSingleString(r))
        r <- summary(r)$description
    cat(class(object), "object\nresource:", r, "\nresolution:", res, "\n")
})
