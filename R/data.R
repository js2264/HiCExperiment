#' @title Example datasets provided in `HiCExperiment` & `HiContactsData`
#' 
#' @name data
#' @aliases centros_yeast
#' @aliases contacts_yeast
#' @aliases contacts_yeast_eco1
#' @aliases full_contacts_yeast
#' 
#' @format An object of class \code{"GRanges"}.
#' @docType data
#' @usage data(centros_yeast)
#' @source HiContacts
#' @examples
#' data(centros_yeast)
#' centros_yeast
#' contacts_yeast()
#' contacts_yeast_eco1()
#' full_contacts_yeast()
NULL 

"centros_yeast"

#' @export

contacts_yeast <- function() {
    fpath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
    import(fpath, 'II', resolution = 1000, format = 'cool')
}

#' @export

contacts_yeast_eco1 <- function() {
    fpath <- HiContactsData::HiContactsData('yeast_eco1', 'mcool')
    import(fpath, 'II', resolution = 1000, format = 'cool')
}

#' @export

full_contacts_yeast <- function() {
    env_ <- new.env(parent = emptyenv())
    data(centros_yeast, envir = env_)
    fpath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
    x <- import(fpath, resolution = 16000, format = 'cool')
    topologicalFeatures(x, 'centromeres') <- env_$centros_yeast
    return(x)
}
