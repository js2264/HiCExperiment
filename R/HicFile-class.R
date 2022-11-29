#' @title `HicFile` S4 class
#' 
#' @name HicFile-class
#' @aliases HicFile
#' 
#' @description
#' 
#' The `HicFile` class describes a `BiocFile` object, pointing to the location 
#' of a .hic file (usually created with juicer) and containing 3 
#' additional slots:
#' 
#' 1. resolution: at which resolution the associated .hic file should be parsed;
#' 2. pairsFile: the path (in plain character) to an optional pairs file 
#'   (stored as a `PairsFile` object);
#' 2. metadata: a list metadata
#' 
#' @param path String; path to a .hic file
#' @param resolution numeric; resolution to use with mcool file
#' @param pairsFile String; path to a pairs file
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
#' hic <- HicFile(hic_path, resolution = 16000, pairsFile = pairs_path)
NULL

#' @export

setClass('HicFile', contains = 'ContactsFile')

#' @export 

HicFile <- function(path, resolution = NULL, pairsFile = NULL, metadata = list()) {
    path <- gsub('~', Sys.getenv('HOME'), path)
    check_hic_file(path)
    if (is_hic(path) & is.null(resolution)) 
        resolution <- lsHicResolutions(path)[1]
    check_hic_format(path, resolution)
    if (!S4Vectors::isSingleString(path))
        stop('"filename" must be a single string, specifiying a path')
    new(
        'HicFile', 
        resource = path, 
        resolution = resolution,
        pairsFile = PairsFile(pairsFile),
        metadata = metadata
    )
}

