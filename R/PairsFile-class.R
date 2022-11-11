#' @title `PairsFile` S4 class and methods
#'
#' @importFrom methods setClass
#' @importClassesFrom BiocIO BiocFile
#' @export
#' @rdname PairsFile 

setClass('PairsFile', contains = 'BiocFile')

#' @rdname PairsFile 
#' @importFrom S4Vectors isSingleString
#' @export
#' @examples 
#' pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#' PairsFile(pairs_path)

PairsFile <- function(path) {
    if (!S4Vectors::isSingleString(path))
        stop('"filename" must be a single string, specifiying a path')
    new('PairsFile', resource = path)
}
