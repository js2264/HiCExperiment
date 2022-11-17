#' @title `HiCoolFile` S4 class
#' 
#' @name HiCoolFile-class
#' @aliases HiCoolFile
#' 
#' @description
#' 
#' The `HiCoolFile` class describes a `BiocFile` object, pointing to the location 
#' of an (m)cool file and containing an additional `resolution` slot, 
#' specifying at which resolution the associated mcool file should be parsed 
#' The `HiCoolFile` class corresponds to an object linking several files
#' generated during the mapping and processing of paired-end sequencing reads
#' by `HiCool::HiCool()`. It extends the `CoolFile` data structure with 
#' two additional slots: 
#' 
#' 1. pairsFile: the path (in plain character) to the pairs file created 
#'   by `HiCool::HiCool()`, stored in a `PairsFile` object;
#' 2. metadata: a list with two elements: `log` (path to `HiCool` processing 
#'   log file) and `stats` (aggregating some stats from `HiCool` mapping).
#'
#' @slot coolFile CoolFile object
#' @slot pairsFile PairsFile object
#' 
#' @param mcool String; path to a (m)cool file
#' @param resolution numeric; resolution to use with mcool file
#' @param pairs String; path to a pairs file
#' @param metadata list; contains path to log file and stats file from hicstuff
#' 
#' @importFrom S4Vectors metadata
#' @importFrom methods setClass
#' @importClassesFrom BiocIO BiocFile
#' @examples
#' mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#' hcf <- HiCoolFile(mcool_path, resolution = 4000, pairs = pairs_path)
NULL

#' @export 

setClass('HiCoolFile', contains = c('CoolFile', 'Annotated'),
    slots = list(
        coolFile = 'CoolFile',
        pairsFile = 'PairsFile'
    )
)

#' @export 

HiCoolFile <- function(mcool, resolution, pairs, metadata = list(log = '', stats = '')) {
    new(
        'HiCoolFile', 
        resource = mcool, 
        resolution = resolution,
        metadata = metadata,
        coolFile = CoolFile(path = mcool, resolution = resolution), 
        pairsFile = PairsFile(pairs)
    )
}
