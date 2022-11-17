#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom InteractionSet GInteractions

setClassUnion("GRangesOrGInteractions", members = c("GRanges", "GInteractions"))
setClassUnion("characterOrNULL", members = c("character", "NULL"))
setClassUnion("numericOrNULL", members = c("numeric", "NULL"))

#' @title `HiCExperiment` S4 class
#' 
#' @name HiCExperiment-class
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
    pairsFile = NULL
) {
    
    ## -- Check that provided file is valid
    if (is(file, 'CoolFile')) {
        resolution <- resolution(file)
        file <- BiocIO::resource(file)
        if (is_mcool(file) & is.null(resolution)) resolution <- lsCoolResolutions(file)[1]
    }
    check_cool_file(file)
    if (is_mcool(file) & is.null(resolution)) resolution <- lsCoolResolutions(file)[1]
    check_cool_format(file, resolution)

    ## -- Read resolutions
    resolutions <- lsCoolResolutions(file)
    if (is_mcool(file)) {
        res <- resolutions[length(resolutions)]
        if (!is.null(resolution)) {
            current_res <- resolution
        }
        else {
            current_res <- res
        }
    } 
    else {
        res <- resolutions[length(resolutions)]
        current_res <- NULL
    }

    ## -- Read interactions
    gis <- cool2gi(file, resolution = current_res, coords = focus)
    mcols <- GenomicRanges::mcols(gis)
    GenomicRanges::mcols(gis) <- mcols[, c('bin_id1', 'bin_id2')]

    ## -- Check pairs file
    if (!is.null(pairsFile)) {
        if (!file.exists(pairsFile)) {
            if (pairsFile == "") {
                pairsFile <- NULL
            }
            else {
                stop("Provided pairsFile does not exist. Aborting now.")
            }
        }
    }

    ## -- Create contact object
    x <- methods::new("HiCExperiment", 
        fileName = as.character(file),
        focus = focus, 
        resolutions = resolutions, 
        resolution = ifelse(is_mcool(file), current_res, res), 
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
