#' @title `ContactsFile` S4 class
#' 
#' @name ContactsFile-class
#' @aliases ContactsFile
#' 
#' @description
#' 
#' The `ContactsFile` class describes a `BiocFile` object, pointing to the location 
#' of an Hi-C matrix file (cool, mcool, hic, hicpro, ...) and containing additional slots:
#' 
#' 1. resolution: at which resolution the associated mcool file should be parsed 
#' 2. pairsFile: the path (in plain character) to an optional pairs file 
#'   (stored as a `PairsFile` object);
#' 3. metadata: a list. If the CoolFile is created by `HiCool`, it will contain 
#'   two elements: `log` (path to `HiCool` processing log file) and `stats` 
#'   (aggregating some stats from `HiCool` mapping).
#'
#' @slot resolution numeric value or NULL 
#' @slot pairsFile PairsFile object
#' @slot metadata list
#' 
#' @param path String; path to an Hi-C matrix file (cool, mcool, hic, hicpro) 
#' @param resolution numeric; resolution to use with Hi-C matrix file
#' @param pairsFile String; path to a pairs file
#' @param metadata list.
#' 
#' @importFrom S4Vectors metadata
#' @importFrom methods setClass
#' @importClassesFrom BiocIO BiocFile
#' @include HiCExperiment-class.R
#' @include PairsFile-class.R

setClass('ContactsFile', contains = c('BiocFile', 'Annotated'),
    slots = list(
        resolution = 'numericOrNULL', 
        pairsFile = 'PairsFileOrNULL'
    )
)
