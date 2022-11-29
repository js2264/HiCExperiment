#' @title `HicproFile` S4 class
#' 
#' @name HicproFile-class
#' @aliases HicproFile
#' 
#' @description
#' 
#' The `HicproFile` class describes a `BiocFile` object, pointing to the location 
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
#' @param path String; path to the HiC-Pro output .matrix file (matrix file)
#' @param bed String; path to the HiC-Pro output .bed file (regions file)
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
#' hicpro_path <- HiContactsData::HiContactsData('yeast_wt', 'hic')
#' pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#' hicpro <- HicproFile(hicpro_path, pairsFile = pairs_path)
NULL

#' @export

setClass('HicproFile', contains = 'ContactsFile', slots = list(
    bed = 'character'
))

#' @export 

HicproFile <- function(path, bed = NULL, pairsFile = NULL, metadata = list()) {
    stopifnot(!is.null(bed) & file.exists(bed))
    check_hicpro_files(path, bed)
    bed1 <- vroom::vroom(
        file = bed, 
        col_names = FALSE, 
        n_max = 1, 
        show_col_types = FALSE, 
        progress = FALSE
    )
    resolution <- (bed1[,3] - bed1[,2])[[1]]
    new(
        'HicproFile', 
        resource = path, 
        bed = bed, 
        resolution = resolution,
        pairsFile = PairsFile(pairsFile),
        metadata = metadata
    )
}

