#' @title `CoolFile` S4 class
#' 
#' @name CoolFile-class
#' @aliases McoolFile-class
#' @aliases CoolFile
#' 
#' @description
#' 
#' The `CoolFile` class describes a `BiocFile` object, pointing to the location 
#' of an (m)cool file and containing additional slots:
#' 
#' 1. resolution: at which resolution the associated mcool file should be parsed 
#' 2. pairsFile: the path (in plain character) to an optional pairs file 
#'   (stored as a `PairsFile` object);
#' 3. metadata: a list with two elements: `log` (path to `HiCool` processing 
#'   log file) and `stats` (aggregating some stats from `HiCool` mapping).
#'
#' @slot resolution numeric value or NULL 
#' @slot pairsFile PairsFile object
#' @slot metadata list
#' 
#' @param path String; path to a (m)cool file
#' @param resolution numeric; resolution to use with mcool file
#' @param pairs String; path to a pairs file
#' @param metadata list; contains path to log file and stats file from hicstuff
#' 
#' @importFrom S4Vectors metadata
#' @importFrom methods setClass
#' @importClassesFrom BiocIO BiocFile
#' @include HiCExperiment-class.R
#' @include PairsFile-class.R
#' 
#' @examples
#' mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#' hcf <- CoolFile(mcool_path, resolution = 4000, pairs = pairs_path)
NULL

#' @export

setClass('CoolFile', contains = c('BiocFile', 'Annotated'),
    slots = list(
        resolution = 'numericOrNULL', 
        pairsFile = 'PairsFile'
    )
)

#' @export

setClass('McoolFile', contains = 'CoolFile')

#' @export 

CoolFile <- function(path, resolution = NULL, pairs = "", metadata = list(log = '', stats = '')) {
    check_cool_file(path)
    if (is_mcool(path) & is.null(resolution)) 
        resolution <- lsCoolResolutions(path)[1]
    check_cool_format(path, resolution)
    if (!S4Vectors::isSingleString(path))
        stop('"filename" must be a single string, specifiying a path')
    new(
        'CoolFile', 
        resource = path, 
        resolution = resolution,
        pairsFile = PairsFile(pairs),
        metadata = metadata
    )
}
