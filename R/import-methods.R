#' @title HiCExperiment import methods
#' 
#' @name import-methods
#' @aliases import
#' @aliases import,CoolFile-method
#' @aliases import,HicFile-method
#' @aliases import,PairsFile-method
#' @aliases import,CoolFile,ANY,ANY-method
#' @aliases import,HicFile,ANY,ANY-method
#' @aliases import,PairsFile,ANY,ANY-method
#' 
#' @description 
#' 
#' Import methods for data structures implemented in the HiCExperiment package. 
#' HiCExperiment package implements methods to faciliate 
#' the import of (m)cool files and pairs files in R, as HiCExperiment  
#' or GenomicInteractions objects.
#' 
#' @param con Path or connection to a cool, mcool or pairs file. Con argument 
#'   can also be a CoolFile or PairsFile object. 
#' @param format The format of the output. If missing and 'con' is a filename, 
#'    the format is derived from the file extension. 
#'    This argument is unnecessary when 'con' is a derivative of 'BiocFile'.
#' @param text If 'con' is missing, this can be a character vector directly 
#'    providing the string data to import.
#' @param ... e.g. `resolution = ...`; parameters to pass to the 
#'    format-specific method.
#' 
#' @usage import(con, format, text, ...)
#' 
#' @importFrom BiocIO import
#' @importFrom BiocGenerics path
#' @exportMethod import 
NULL

#' @export

setMethod('import', 'CoolFile', function(con, ...) {

    path <- BiocGenerics::path(con)
    path <- gsub('~', Sys.getenv('HOME'), path)
    stopifnot(file.exists(path))
    if (!check_cool_file(path)) stop("Provided file is not a valid .(m)cool file")
    params <- list(...)

    ## -- Check parsed arguments
    if ('resolution' %in% names(params)) {
        check_cool_format(path, params[['resolution']])
        resolution <- params[['resolution']]
    } else {
        if (is_mcool(path)) {
            resolution <- lsCoolResolutions(path)[1]
        }
        else {
            resolution <- NULL
        }
    }
    if ('focus' %in% names(params)) {
        focus <- params[['focus']]
    } else {
        focus <- NULL
    }
    if ('pairsFile' %in% names(params)) {
        if (!file.exists(params[['pairsFile']])) {
            if (params[['pairsFile']] == "") {
                pairsFile <- NULL
            }
            else {
                stop("Provided pairsFile does not exist. Aborting now.")
            }
        } else {
            pairsFile <- params[['pairsFile']]
        }
    } else {
        pairsFile <- NULL
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
        metadata <- list()
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
    params <- list(...)

    ## -- Check parsed arguments
    if ('resolution' %in% names(params)) {
        # check_hic_format(path, params[['resolution']])
        resolution <- params[['resolution']]
    } else {
        resolution <- strawr::readHicBpResolutions(path)[1]
    }
    if ('focus' %in% names(params)) {
        focus <- params[['focus']]
        focus <- gsub("-", ":", params[['focus']])
    } else {
        stop("`focus` argument is required for .hic files.")
    }
    if ('pairsFile' %in% names(params)) {
        if (!file.exists(params[['pairsFile']])) {
            if (params[['pairsFile']] == "") {
                pairsFile <- NULL
            }
            else {
                stop("Provided pairsFile does not exist. Aborting now.")
            }
        } else {
            pairsFile <- params[['pairsFile']]
        }
    } else {
        pairsFile <- NULL
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
        metadata <- list()
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

setMethod('import', signature = 'PairsFile', function(con, ...) {

    con <- BiocGenerics::path(con)
    stopifnot(file.exists(con))
    pairs2gi(con, ...)

})
