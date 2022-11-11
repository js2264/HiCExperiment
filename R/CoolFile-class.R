setClassUnion("numericOrNULL", members = c("numeric", "NULL"))

#' @title `CoolFile` S4 class and methods
#'
#' @slot resolution numeric value or NULL 
#' 
#' @importFrom methods setClass
#' @importClassesFrom BiocIO BiocFile
#' @export
#' @rdname CoolFile 

setClass('CoolFile', contains = 'BiocFile',
    slots = list(resolution = 'numericOrNULL')
)

#' @rdname CoolFile 
#' @importFrom S4Vectors isSingleString
#' @export
#' @examples 
#' cool_path <- HiContactsData::HiContactsData('yeast_wt', 'cool')
#' CoolFile(cool_path)
#' cool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' CoolFile(cool_path, resolution = 1000)

CoolFile <- function(path, resolution = NULL) {
    check_cool_file(path)
    check_cool_format(path, resolution)
    if (!S4Vectors::isSingleString(path))
        stop('"filename" must be a single string, specifiying a path')
    new('CoolFile', resource = path, resolution = resolution)
}
