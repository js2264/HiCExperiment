#' @title Example datasets provided in `HiCExperiment` & `HiContactsData`
#' 
#' @name data
#' @aliases centros_yeast
#' @aliases contacts_yeast
#' @aliases contacts_yeast_eco1
#' 
#' @format An object of class \code{"GRanges"}.
#' @docType data
#' @usage data(centros_yeast)
#' @source HiContacts
#' @param full Whether to import all interactions
#' @examples
#' data(centros_yeast)
#' centros_yeast
#' contacts_yeast()
NULL 

"centros_yeast"

#' @export
#' @rdname data

contacts_yeast <- function(full = FALSE) {
    env_ <- new.env(parent = emptyenv())
    data(centros_yeast, envir = env_)
    fpath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
    if (full) {
        x <- import(fpath, resolution = 16000, format = 'cool')
    }
    else {
        x <- import(fpath, 'II', resolution = 16000, format = 'cool')
    }
    topologicalFeatures(x, 'centromeres') <- env_$centros_yeast
    return(x)
}

#' @export
#' @rdname data

contacts_yeast_eco1 <- function(full = FALSE) {
    fpath <- HiContactsData::HiContactsData('yeast_eco1', 'mcool')
    import(fpath, 'II', resolution = 16000, format = 'cool')
}
