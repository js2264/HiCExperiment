setClassUnion("numericOrNULL", members = c("numeric", "NULL"))

#' @title `CoolFile` S4 class
#' 
#' @name CoolFile-class
#' @aliases CoolFile
#' 
#' @description
#' 
#' The `CoolFile` class describes a `BiocFile` object, pointing to the location 
#' of an (m)cool file and containing an additional `resolution` slot, 
#' specifying at which resolution the associated mcool file should be parsed 
#' when imported in R. 
#'
#' @slot resolution numeric value or NULL 
#' 
#' @param path String; path to a (m)cool file
#' @param resolution numeric; resolution to use with an mcool file
#' 
#' @importFrom methods setClass
#' @importClassesFrom BiocIO BiocFile
#' @importFrom S4Vectors isSingleString
#' @examples 
#' cool_path <- HiContactsData::HiContactsData('yeast_wt', 'cool')
#' CoolFile(cool_path)
#' cool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' cf <- CoolFile(cool_path, resolution = 4000)
NULL

#' @export

setClass('CoolFile', contains = 'BiocFile',
    slots = list(resolution = 'numericOrNULL')
)

#' @export

setClass('McoolFile', contains = 'CoolFile')

#' @export

CoolFile <- function(path, resolution = NULL) {
    check_cool_file(path)
    if (is_mcool(path) & is.null(resolution)) 
        resolution <- lsCoolResolutions(path)[1]
    check_cool_format(path, resolution)
    if (!S4Vectors::isSingleString(path))
        stop('"filename" must be a single string, specifiying a path')
    new('CoolFile', resource = path, resolution = resolution)
}
