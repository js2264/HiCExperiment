#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom InteractionSet GInteractions

setClassUnion("GRangesOrGInteractions", members = c("GRanges", "GInteractions"))
setClassUnion("characterOrNULL", members = c("character", "NULL"))
setClassUnion("numericOrNULL", members = c("numeric", "NULL"))

#' @title `HiCExperiment` S4 class
#' 
#' @name HiCExperiment
#' @rdname HiCExperiment
#' 
#' @description
#' 
#' The `HiCExperiment` class describes (m)cool files imported in R, either 
#' through the `HiCExperiment` constructor function or using the `import` 
#' method implemented by `HiCExperiment` package. 
#'
#' @slot fileName Path of (m)cool file
#' @slot focus Chr. coordinates for which interaction counts are extracted 
#'   from the .(m)cool file.
#' @slot resolutions Resolutions available in the .(m)cool file.
#' @slot resolution Current resolution
#' @slot interactions Genomic Interactions extracted from the .(m)cool object
#' @slot scores Available interaction scores. 
#' @slot topologicalFeatures Topological features associated with the dataset 
#'   (e.g. loops (\<Pairs\>), borders (\<GRanges\>), 
#'   viewpoints (\<GRanges\>), etc...)
#' @slot pairsFile Path to the .pairs file associated with the .(m)cool file
#' @slot metadata metadata associated with the .(m)cool file.
#' 
#' @param file CoolFile or plain path to a (m)cool file
#' @param resolution Resolution to use with mcool file
#' @param focus Chromosome coordinates for which 
#'   interaction counts are extracted from the .(m)cool file, provided
#'   as a character string (e.g. "II:4000-5000"). If not provided, 
#'   the entire (m)cool file will be imported. 
#' @param metadata list of metadata
#' @param topologicalFeatures topologicalFeatures provided as a named SimpleList
#' @param pairsFile Path to an associated .pairs file
#' @param ... Extra arguments
#' 
#' @return An `HiCExperiment` object.
#' 
#' @importFrom methods setClass
#' @importClassesFrom S4Vectors Annotated
#' @examples 
#' library(HiCExperiment)
#' mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' HiCExperiment(mcool_path, resolution = 16000)
NULL

#' @rdname HiCExperiment
#' @export

methods::setClass("HiCExperiment", 
    contains = c("Annotated"), 
    slots = c(
        fileName = "character",
        focus = "characterOrNULL", 
        resolutions = "numeric", 
        resolution = "numeric", 
        interactions = "GInteractions",
        scores = "SimpleList", 
        topologicalFeatures = "SimpleList",
        pairsFile = "characterOrNULL",
        metadata = "list"
    )
)

#' @rdname HiCExperiment
#' @export

 HiCExperiment <- function(
    file, 
    resolution = NULL, 
    focus = NULL, 
    metadata = list(), 
    topologicalFeatures = S4Vectors::SimpleList(
        'loops' = S4Vectors::Pairs(
            GenomicRanges::GRanges(), 
            GenomicRanges::GRanges()
        ), 
        'borders' = GenomicRanges::GRanges(), 
        'compartments' = GenomicRanges::GRanges(), 
        'viewpoints' = GenomicRanges::GRanges()
    ), 
    pairsFile = NULL, 
    ...
) {
    file <- gsub('~', Sys.getenv('HOME'), file)
    stopifnot(file.exists(file))
    params <- list(...)
    if ("bed" %in% names(params)) bed <- params[['bed']]

    if (is_cool(file) | is_mcool(file)) {
        return(.HiCExperimentFromCoolFile(
            file = file,
            resolution = resolution,
            focus = focus,
            metadata = metadata,
            topologicalFeatures = topologicalFeatures,
            pairsFile = pairsFile
        ))
    }
    if (is_hic(file)) {
        return(.HiCExperimentFromHicFile(
            file = file,
            resolution = resolution,
            focus = focus,
            metadata = metadata,
            topologicalFeatures = topologicalFeatures,
            pairsFile = pairsFile
        ))
    }
    if (is_hicpro_matrix(file) & is_hicpro_regions(bed)) {
        return(.HiCExperimentFromHicproFile(
            file = file,
            bed = bed,
            metadata = metadata,
            topologicalFeatures = topologicalFeatures,
            pairsFile = pairsFile
        ))
    }
    return(NA)
}

setValidity("HiCExperiment",
    function(object) {
        if (!is(focus(object), "characterOrNULL"))
            return("'focus' slot must be a characterOrNULL")
        if (!is(resolutions(object), "numeric"))
            return("'resolutions' slot must be an numeric vector")
        if (!is(scores(object), "SimpleList"))
            return("'scores' slot must be a SimpleList")
        TRUE
    }
)

#' @rdname HiCExperiment

.HiCExperimentFromCoolFile <- function(
    file, 
    resolution = NULL, 
    focus = NULL, 
    metadata = list(), 
    topologicalFeatures = S4Vectors::SimpleList(
        'loops' = S4Vectors::Pairs(
            GenomicRanges::GRanges(), 
            GenomicRanges::GRanges()
        ), 
        'borders' = GenomicRanges::GRanges(), 
        'compartments' = GenomicRanges::GRanges(), 
        'viewpoints' = GenomicRanges::GRanges()
    ), 
    pairsFile = NULL
) {
    
    ## -- Check that provided file is valid
    file <- gsub('~', Sys.getenv('HOME'), file)
    check_cool_file(file)
    check_cool_format(file, resolution)

    ## -- Read interactions
    gis <- .cool2gi(file, resolution = resolution, coords = focus) |> sort()
    mcols <- GenomicRanges::mcols(gis)
    GenomicRanges::mcols(gis) <- mcols[, c('bin_id1', 'bin_id2')]

    ## -- Create contact object
    x <- methods::new("HiCExperiment", 
        fileName = as.character(file),
        focus = focus, 
        resolutions = lsCoolResolutions(file), 
        resolution = ifelse(is.null(resolution), lsCoolResolutions(file)[1], resolution), 
        interactions = gis, 
        scores = S4Vectors::SimpleList(
            'raw' = as.numeric(mcols$count),
            'balanced' = as.numeric(mcols$score)
        ), 
        topologicalFeatures = topologicalFeatures, 
        pairsFile = pairsFile, 
        metadata = metadata
    )
    methods::validObject(x)
    return(x)
} 

#' @rdname HiCExperiment

.HiCExperimentFromHicFile <- function(
    file, 
    resolution = NULL, 
    focus = NULL, 
    metadata = list(), 
    topologicalFeatures = S4Vectors::SimpleList(
        'loops' = S4Vectors::Pairs(
            GenomicRanges::GRanges(), 
            GenomicRanges::GRanges()
        ), 
        'borders' = GenomicRanges::GRanges(), 
        'compartments' = GenomicRanges::GRanges(), 
        'viewpoints' = GenomicRanges::GRanges()
    ), 
    pairsFile = NULL
) {
    
    ## -- Check that provided file is valid
    file <- gsub('~', Sys.getenv('HOME'), file)
    check_hic_file(file)
    check_hic_format(file, resolution)

    ## -- Read interactions
    gis <- .hic2gi(file, resolution = resolution, coords = focus) |> sort()
    mcols <- GenomicRanges::mcols(gis)
    GenomicRanges::mcols(gis) <- mcols[, c('bin_id1', 'bin_id2')]

    ## -- Create HiCExperiment
    if (!is.null(focus)) focus <- gsub("(:.*[^:]*):", "\\1-", focus)
    x <- methods::new("HiCExperiment", 
        fileName = as.character(file),
        focus = focus, 
        resolutions = rev(strawr::readHicBpResolutions(file)), 
        resolution = resolution, 
        interactions = gis, 
        scores = S4Vectors::SimpleList(
            'raw' = as.numeric(mcols$count),
            'balanced' = as.numeric(mcols$score)
        ), 
        topologicalFeatures = topologicalFeatures, 
        pairsFile = pairsFile, 
        metadata = metadata
    )
    methods::validObject(x)
    return(x)
} 

#' @rdname HiCExperiment

.HiCExperimentFromHicproFile <- function(
    file, 
    bed, 
    metadata = list(), 
    topologicalFeatures = S4Vectors::SimpleList(
        'loops' = S4Vectors::Pairs(
            GenomicRanges::GRanges(), 
            GenomicRanges::GRanges()
        ), 
        'borders' = GenomicRanges::GRanges(), 
        'compartments' = GenomicRanges::GRanges(), 
        'viewpoints' = GenomicRanges::GRanges()
    ), 
    pairsFile = NULL
) {
    
    ## -- Check that provided file is valid
    file <- gsub('~', Sys.getenv('HOME'), file)
    check_hicpro_files(file, bed)

    ## -- Read interactions
    gis <- .hicpro2gi(file, bed) |> sort()
    mcols <- GenomicRanges::mcols(gis)
    GenomicRanges::mcols(gis) <- mcols[, c('bin_id1', 'bin_id2')]
    res <- GenomicRanges::width(regions(gis))[[1]]

    ## -- Create HiCExperiment
    x <- methods::new("HiCExperiment", 
        fileName = as.character(file),
        focus = NULL, 
        resolutions = res, 
        resolution = res, 
        interactions = gis, 
        scores = S4Vectors::SimpleList(
            'counts' = as.numeric(mcols$count)
        ), 
        topologicalFeatures = topologicalFeatures, 
        pairsFile = pairsFile, 
        metadata = c(list(regions = bed), metadata)
    )
    methods::validObject(x)
    return(x)
} 




