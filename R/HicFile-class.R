#' @title `HicFile` S4 class
#' 
#' @name HicFile-class
#' @aliases HicFile
#' 
#' @description
#' 
#' The `HicFile` class describes a `BiocFile` object, pointing to the location 
#' of a .hic file and containing 3 additional slots:
#' 
#' 1. resolution: at which resolution the associated mcool file should be parsed 
#' 2. pairsFile: the path (in plain character) to an optional pairs file 
#'   (stored as a `PairsFile` object);
#' 2. metadata: a list metadata
#'
#' @slot resolution numeric value or NULL 
#' @slot metadata list
#' 
#' @param path String; path to a .hic file
#' @param resolution numeric; resolution to use with mcool file
#' @param pairs String; path to a pairs file
#' @param metadata list.
#' 
#' @importFrom S4Vectors metadata
#' @importFrom methods setClass
#' @importClassesFrom BiocIO BiocFile
#' @include HiCExperiment-class.R
#' @include PairsFile-class.R
#' 
#' @examples
#' hic_path <- HiContactsData::HiContactsData('yeast_wt', 'hic')
#' pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#' hic <- HicFile(hic_path, resolution = pairs = pairs_path)
NULL

#' @export

setClass('HicFile', contains = c('BiocFile', 'Annotated'),
    slots = list(
        resolution = 'numericOrNULL', 
        pairsFile = 'PairsFile'
    )
)

#' @export 

HicFile <- function(path, resolution = NULL, pairs = "", metadata = list()) {
    path <- gsub('~', Sys.getenv('HOME'), path)
    # check_hic_file(path)
    if (is.null(resolution)) 
        resolution <- strawr::readHicBpResolutions(path)[1]
    if (!S4Vectors::isSingleString(path))
        stop('"filename" must be a single string, specifiying a path')
    new(
        'HicFile', 
        resource = path, 
        resolution = resolution,
        pairsFile = PairsFile(pairs),
        metadata = metadata
    )
}

