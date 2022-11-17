#' @title Example datasets provided in `HiCExperiment` & `HiContactsData`
#' 
#' @name data
#' 
#' @format An object of class \code{"GRanges"}.
#' @docType data
#' @usage data(centros_yeast)
#' @source HiContacts
#' @importFrom HiContactsData HiContactsData
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
    HiCExperiment(fpath, 'II', resolution = 1000)
}

#' @export

contacts_yeast_eco1 <- function() {
    fpath <- HiContactsData::HiContactsData('yeast_eco1', 'mcool')
    HiCExperiment(fpath, 'II', resolution = 1000)
}

#' @export

full_contacts_yeast <- function() {
    env_ <- new.env(parent = emptyenv())
    data(centros_yeast, envir = env_)
    fpath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
    x <- HiCExperiment(fpath, resolution = 16000)
    topologicalFeatures(x, 'centromeres') <- env_$centros_yeast
    return(x)
}
