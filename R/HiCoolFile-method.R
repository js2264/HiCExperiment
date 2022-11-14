#' @rdname HiCoolFile
#'
#' @name show
#' @docType methods
#' @aliases show,HiCoolFile-method
#'
#' @param object A \code{HiCoolFile} object.
#'
#' @export
#' @examples 
#' show(contacts_yeast)

setMethod("show", signature("HiCoolFile"), function(object) {
    r <- BiocIO::resource(object)
    res <- resolution(object)
    if (is.null(res)) res = lsCoolResolutions(r)[1]
    if (!S4Vectors::isSingleString(r))
        r <- summary(r)$description
    cat(
        class(object), "object\nresource:", r, 
        "\nresolution:", res, 
        "\npairsFile:", pairsFile(object), "\n"
        )
})
