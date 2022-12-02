#' @title HiCExperiment import methods
#' 
#' @name import-methods
#' @aliases import
#' @aliases import,CoolFile-method
#' @aliases import,HicFile-method
#' @aliases import,HicproFile-method
#' @aliases import,PairsFile-method
#' @aliases import,CoolFile,ANY,ANY-method
#' @aliases import,HicFile,ANY,ANY-method
#' @aliases import,HicproFile,ANY,ANY-method
#' @aliases import,PairsFile,ANY,ANY-method
#' 
#' @description 
#' 
#' Import methods for data structures implemented in the HiCExperiment package. 
#' HiCExperiment package implements methods to faciliate 
#' the import of Hi-C files (.(m)cool, .hic, HiC-Pro derived) 
#' and pairs files in R, as HiCExperiment or GenomicInteractions objects.
#' 
#' @param con Path or connection to a cool, mcool, .hic or HiC-Pro derived files. 
#'   Can also be a path to a pairs file. 
#' @param format The format of the output. If missing and 'con' is a filename, 
#'    the format is derived from the file extension. 
#'    This argument is unnecessary when files are directly provided as 
#'    `CoolFile`, `HicFile`, `HicproFile` or `PairsFile`.
#' @param text If 'con' is missing, this can be a character vector directly 
#'    providing the string data to import.
#' @param ... e.g. `resolution = ...`; parameters to pass to  
#'    format-specific methods.
#' 
#' @usage import(con, format, text, ...)
#' 
#' @importFrom BiocIO import
#' @importFrom BiocGenerics path
#' @exportMethod import 
#' @examples
#' # ---- Importing .(m)cool files 
#' mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' import(mcool_path, resolution = 16000, format = 'cool')
#' 
#' # ---- Importing .hic files 
#' hic_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' import(hic_path, resolution = 2000, focus = 'II', format = 'hic')
#' 
#' # ---- Importing HiC-Pro files 
#' hicpro_matrix_path <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_matrix')
#' hicpro_bed_path <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_bed')
#' import(hicpro_matrix_path, bed = hicpro_bed_path, format = 'hicpro')
NULL

#' @export

setMethod('import', 'CoolFile', function(con, ...) {

    path <- BiocGenerics::path(con)
    path <- gsub('~', Sys.getenv('HOME'), path)
    stopifnot(file.exists(path))
    if (!check_cool_file(path)) stop("Provided file is not a valid .(m)cool file")

    ## -- Handle parsed arguments. Priority is given to 
    #       explicitly provided arguments, then to resolution/pairsFile/metadata 
    #       stored in *File.
    params <- list(...)
    if ('resolution' %in% names(params)) {
        check_cool_format(path, params[['resolution']])
        resolution <- as.integer(params[['resolution']])
    } else {
        resolution <- resolution(con)
    } 
    if ('focus' %in% names(params)) {
        focus <- params[['focus']]
    } else {
        focus <- NULL
    }
    if ('pairsFile' %in% names(params)) {
        if (!file.exists(params[['pairsFile']])) {
            stop("Provided pairsFile does not exist. Aborting now.")
        } else {
            pairsFile <- params[['pairsFile']]
        }
    } else {
        pairsFile <- pairsFile(con)
    }
    if ('topologicalFeatures' %in% names(params)) {
        topologicalFeatures <- params[['topologicalFeatures']]
    } else {
        topologicalFeatures <- S4Vectors::SimpleList(
            loops = S4Vectors::Pairs(GenomicRanges::GRanges(), GenomicRanges::GRanges()), 
            borders = GenomicRanges::GRanges(), 
            compartments = GenomicRanges::GRanges(), 
            viewpoints = GenomicRanges::GRanges()
        )
    }
    if ('metadata' %in% names(params)) {
        metadata <- params[['metadata']]
    } else {
        metadata <- metadata(con)
    }

    ## -- Create HiCExperiment
    .HiCExperimentFromCoolFile(
        path,
        resolution = resolution,
        focus = focus,
        metadata = metadata,
        topologicalFeatures = topologicalFeatures,
        pairsFile = pairsFile
    )

})

#' @export

setMethod('import', 'HicFile', function(con, ...) {

    path <- BiocGenerics::path(con)
    path <- gsub('~', Sys.getenv('HOME'), path)
    stopifnot(file.exists(path))
    if (!check_hic_file(path)) stop("Provided file is not a valid .hic file")

    ## -- Handle parsed arguments. Priority is given to 
    #       explicitly provided arguments, then to resolution/pairsFile/metadata 
    #       potentially stored in *File.
    params <- list(...)
    if ('resolution' %in% names(params)) {
        check_hic_format(path, params[['resolution']])
        resolution <- as.integer(params[['resolution']])
    } else {
        resolution <- resolution(con)
    }
    if ('focus' %in% names(params)) {
        focus <- params[['focus']]
        focus <- gsub("-", ":", params[['focus']])
    } else {
        focus <- NULL
    }
    if ('pairsFile' %in% names(params)) {
        if (!file.exists(params[['pairsFile']])) {
            stop("Provided pairsFile does not exist. Aborting now.")
        } else {
            pairsFile <- params[['pairsFile']]
        }
    } else {
        pairsFile <- pairsFile(con)
    }
    if ('topologicalFeatures' %in% names(params)) {
        topologicalFeatures <- params[['topologicalFeatures']]
    } else {
        topologicalFeatures <- S4Vectors::SimpleList(
            loops = S4Vectors::Pairs(GenomicRanges::GRanges(), GenomicRanges::GRanges()), 
            borders = GenomicRanges::GRanges(), 
            compartments = GenomicRanges::GRanges(), 
            viewpoints = GenomicRanges::GRanges()
        )
    }
    if ('metadata' %in% names(params)) {
        metadata <- params[['metadata']]
    } else {
        metadata <- metadata(con)
    }

    ## -- Create HiCExperiment
    .HiCExperimentFromHicFile(
        path,
        resolution = resolution,
        focus = focus,
        metadata = metadata,
        topologicalFeatures = topologicalFeatures,
        pairsFile = pairsFile
    )

})

#' @export

setMethod('import', 'HicproFile', function(con, ...) {

    params <- list(...)
    path <- BiocGenerics::path(con)
    path <- gsub('~', Sys.getenv('HOME'), path)
    stopifnot(file.exists(path))

    ## -- Handle parsed arguments. Priority is given to 
    #       explicitly provided arguments, then to resolution/pairsFile/metadata 
    #       potentially stored in *File.
    if ('bed' %in% names(params)) {
        if (!file.exists(params[['bed']])) {
            stop("Provided regions file does not exist.")
        } else {
            bed <- params[['bed']]
        }
    } else {
        bed <- con@bed
    }
    if (is.null(bed)) stop("Regions file not provided. ")
    if ('pairsFile' %in% names(params)) {
        if (!file.exists(params[['pairsFile']])) {
            stop("Provided pairsFile does not exist.")
        } else {
            pairsFile <- params[['pairsFile']]
        }
    } else {
        pairsFile <- pairsFile(con)
    }
    if ('topologicalFeatures' %in% names(params)) {
        topologicalFeatures <- params[['topologicalFeatures']]
    } else {
        topologicalFeatures <- S4Vectors::SimpleList(
            loops = S4Vectors::Pairs(GenomicRanges::GRanges(), GenomicRanges::GRanges()), 
            borders = GenomicRanges::GRanges(), 
            compartments = GenomicRanges::GRanges(), 
            viewpoints = GenomicRanges::GRanges()
        )
    }
    if ('metadata' %in% names(params)) {
        metadata <- params[['metadata']]
    } else {
        metadata <- metadata(con)
    }

    if (!check_hicpro_files(path, bed))
        stop("Provided file is not a valid HiC-Pro derived file")

    ## -- Create HiCExperiment
    .HiCExperimentFromHicproFile(
        path,
        bed = bed,
        metadata = metadata,
        topologicalFeatures = topologicalFeatures,
        pairsFile = pairsFile
    )

})

#' @export

setMethod('import', 'PairsFile', function(con, ...) {

    con <- BiocGenerics::path(con)
    stopifnot(file.exists(con))
    pairs2gi(con, ...)

})
