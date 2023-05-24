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
#' @aliases availableResolutions,ANY-method
#' @aliases availableResolutions,CoolFile-method
#' @aliases availableResolutions,HicFile-method
#' @aliases availableResolutions,HicproFile-method
#' @aliases availableChromosomes,ANY-method
#' @aliases availableChromosomes,CoolFile-method
#' @aliases availableChromosomes,HicFile-method
#' @aliases availableChromosomes,HicproFile-method
#' 
#' @description 
#' 
#' Import methods to parse Hi-C files (`.(m)cool`, `.hic`, HiC-Pro derived 
#' matrices, pairs files) into data structures implemented in the 
#' HiCExperiment package. 
#' 
#' @param con,x Path or connection to a cool, mcool, .hic or HiC-Pro derived files. 
#'   Can also be a path to a pairs file. 
#' @param format The format of the output. If missing and 'con' is a filename, 
#'    the format is derived from the file extension. 
#'    This argument is unnecessary when files are directly provided as 
#'    `CoolFile`, `HicFile`, `HicproFile` or `PairsFile`.
#' @param text If 'con' is missing, this can be a character vector directly 
#'    providing the string data to import.
#' @param ... Extra parameters to pass to  
#'    format-specific methods. A list of possible arguments is 
#'    provided in the next section.
#' 
#' @section import arguments for ContactFile class:
#' `ContactFile` class gathers `CoolFile`, `HicFile` and `HicproFile` classes. 
#' When importing a `ContactFile` object in R, two main arguments can be 
#' provided besides the `ContactFile` itself: 
#' 
#' - `resolution`:
#'    Resolutions available in the disk-stored contact matrix can be 
#'    listed using `availableResolutions(file)`
#' - `focus`:
#'    A genomic locus (or pair of loci) provided as a string. It can be any 
#'    of the following string structures: 
#' 
#'      - `"II"` or `"II:20001-30000"`: this will extract a symmetrical 
#'      square HiCExperiment object, of an entire chromosome or an portion of it.
#'      - `"II|III"` or `"II:20001-30000|III:40001-90000"`: 
#'      this will extract a non-symmetrical HiCExperiment object, 
#'      with an entire or portion of different chromosomes on each axis. 
#' 
#' @usage import(con, format, text, ...)
#' @return A `HiCExperiment` or `GInteractions` object
#' 
#' @importFrom BiocIO import
#' @importFrom BiocGenerics path
#' @examples
#' ################################################################
#' ## ----------- Importing .(m)cool contact matrices ---------- ##
#' ################################################################
#' 
#' mcoolPath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' availableResolutions(mcoolPath)
#' availableChromosomes(mcoolPath)
#' import(mcoolPath, resolution = 16000, focus = 'XVI', format = 'cool')
#' 
#' ################################################################
#' ## ------------ Importing .hic contact matrices ------------- ##
#' ################################################################
#' 
#' hicPath <- HiContactsData::HiContactsData('yeast_wt', 'hic')
#' availableResolutions(hicPath)
#' availableChromosomes(hicPath)
#' import(hicPath, resolution = 16000, focus = 'XVI', format = 'hic')
#' 
#' ################################################################
#' ## ------- Importing HiC-Pro derived contact matrices ------- ##
#' ################################################################
#' 
#' hicproMatrixPath <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_matrix')
#' hicproBedPath <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_bed')
#' availableResolutions(hicproMatrixPath, hicproBedPath)
#' availableChromosomes(hicproMatrixPath, hicproBedPath)
#' import(hicproMatrixPath, bed = hicproBedPath, format = 'hicpro')
NULL

#' @exportMethod import 

setMethod('import', 'CoolFile', function(con, ...) {

    path <- BiocGenerics::path(con)
    path <- gsub('~', Sys.getenv('HOME'), path)
    stopifnot(file.exists(path))
    if (!.check_cool_file(path)) stop("Provided file is not a valid .(m)cool file")

    ## -- Handle parsed arguments. Priority is given to 
    #       explicitly provided arguments, then to resolution/pairsFile/metadata 
    #       stored in *File.
    params <- list(...)
    if ('resolution' %in% names(params)) {
        .check_cool_format(path, params[['resolution']])
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
            compartments = GenomicRanges::GRanges(), 
            borders = GenomicRanges::GRanges(), 
            loops = InteractionSet::GInteractions(GenomicRanges::GRanges(), GenomicRanges::GRanges()), 
            viewpoints = GenomicRanges::GRanges()
        )
    }
    if ('metadata' %in% names(params)) {
        metadata <- params[['metadata']]
    } else {
        metadata <- metadata(con)
    }
    if ('compute.detrend' %in% names(params)) {
        if (params[['detrend']]) {
            metadata[['detrending_model']] <- detrendingModel(path, resolution)
        }
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
    if (!.check_hic_file(path)) stop("Provided file is not a valid .hic file")

    ## -- Handle parsed arguments. Priority is given to 
    #       explicitly provided arguments, then to resolution/pairsFile/metadata 
    #       potentially stored in *File.
    params <- list(...)
    if ('resolution' %in% names(params)) {
        .check_hic_format(path, params[['resolution']])
        resolution <- as.integer(params[['resolution']])
    } else {
        resolution <- resolution(con)
    }
    if ('focus' %in% names(params)) {
        focus <- params[['focus']]
        # focus <- gsub("-", ":", params[['focus']])
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
            compartments = GenomicRanges::GRanges(), 
            borders = GenomicRanges::GRanges(), 
            loops = InteractionSet::GInteractions(GenomicRanges::GRanges(), GenomicRanges::GRanges()), 
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
            stop("Provided regions file (`bed` argument) does not exist.")
        } else {
            bed <- params[['bed']]
        }
    } else {
        bed <- con@bed
    }
    if (is.null(bed)) stop("Regions file (`bed` argument) not provided. ")
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
            compartments = GenomicRanges::GRanges(), 
            borders = GenomicRanges::GRanges(), 
            loops = InteractionSet::GInteractions(GenomicRanges::GRanges(), GenomicRanges::GRanges()), 
            viewpoints = GenomicRanges::GRanges()
        )
    }
    if ('metadata' %in% names(params)) {
        metadata <- params[['metadata']]
    } else {
        metadata <- metadata(con)
    }
    if ('focus' %in% names(params)) {
        warning('HiC-Pro contact matrix does not support random access. `focus` argument is ignored.')
    }

    if (!.check_hicpro_files(path, bed))
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
    .pairs2gi(con, ...)

})
