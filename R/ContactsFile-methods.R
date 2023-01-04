#' @title `ContactsFile` methods
#' 
#' @name ContactsFile-methods
#' @rdname ContactsFile-class
#' @aliases pairsFile,ContactsFile-method
#' @aliases resolution,ContactsFile-method
#' @aliases metadata<-,ContactsFile-method
#' @aliases metadata<-,ContactsFile,list-method
#' 
#' @description
#' 
#' ContactsFile methods.
#'
#' @param object A \code{ContactsFile} object.
#' @param x A \code{ContactsFile} object.
#'
#' @importFrom BiocGenerics path
#' @include ContactsFile-class.R 
#' @include PairsFile-class.R 
NULL

#' @export

setMethod("pairsFile", "ContactsFile", function(x) {
    if (is.null(x@pairsFile)) return(NULL)
    BiocGenerics::path(x@pairsFile)
})

#' @export

setMethod("resolution", "ContactsFile", function(x) x@resolution)

#' @export

setMethod("metadata<-", signature(x = "ContactsFile", value = "list"), function(x, value) {
    x@metadata <- value
    x
})
