#' @title `AggrHiCExperiment` S4 class
#' 
#' @name AggrHiCExperiment
#' @rdname AggrHiCExperiment
#' 
#' @description
#' 
#' The `AggrHiCExperiment` extends `HiCExperiment` class
#'
#' @slot fileName Path of Hi-C contact file
#' @slot resolutions Resolutions available in the Hi-C contact file.
#' @slot resolution Current resolution
#' @slot interactions Genomic Interactions extracted from the Hi-C contact file
#' @slot scores Available interaction scores. 
#' @slot topologicalFeatures Topological features associated with the dataset 
#'   (e.g. loops (\<Pairs\>), borders (\<GRanges\>), 
#'   viewpoints (\<GRanges\>), etc...)
#' @slot pairsFile Path to the .pairs file associated with the Hi-C contact file
#' @slot metadata metadata associated with the Hi-C contact file.
#' 
#' @include HiCExperiment-class.R
#' 
NULL

#' @rdname AggrHiCExperiment
#' @export

methods::setClass("AggrHiCExperiment", contains = c("HiCExperiment"))

#' @rdname AggrHiCExperiment
#' @export
#' @examples
#' library(rtracklayer)
#' devtools::load_all('.')
#' file = "/home/rsg/.cache/R/fourDNData/78685597e138_4DNFI9FVHJZQ.mcool"
#' file = "/data/20221214_HiContacts_compartments-TADs-loops/data/test.mcool"
#' resolution = 10000 
#' BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)
#' coords_list = import('/home/rsg/Projects/20221214_HiContacts_compartments-TADs-loops/data/test_CTCF.bed.gz')
#' CTCF_bw <- import('/home/rsg/Projects/20221214_HiContacts_compartments-TADs-loops/data/test_CTCF.bigWig', as = 'NumericList', selection = BigWigSelection(resize(coords_list, fix = 'center', width = 250)))
#' coords_list = import('/home/rsg/Projects/20221214_HiContacts_compartments-TADs-loops/data/test_CTCF.bed.gz')
#' coords_list$scoreChIP <- mean(CTCF_bw)
#' snippets <- coords_list[coords_list$score >= quantile(coords_list$score , 0.75) & coords_list$scoreChIP >= quantile(coords_list$scoreChIP , 0.75)]
#' snippets <- GenomicRanges::resize(snippets, width = resolution*100, fix = 'center')
#' snippets <- snippets[order(snippets$scoreChIP, decreasing = TRUE)]
#' pairs <- Pairs(
#'     rep(snippets, each = length(snippets)),
#'     rep(snippets, length(snippets))
#' )
#' d <- distance(pairs)
#' pairs <- pairs[d > 200000 & d < 1000000 & !is.na(d)]

AggrHiCExperiment <- function(
    file, 
    resolution = NULL, 
    snippets, 
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
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
    ...
) {
    file <- gsub('~', Sys.getenv('HOME'), file)
    stopifnot(file.exists(file))
    params <- list(...)
    bed <- NULL ; if ("bed" %in% names(params)) bed <- params[['bed']]
    if (!is.null(resolution)) resolution <- as.integer(resolution)
    
    if (is_cool(file) | is_mcool(file)) {
        return(.AggrHiCExperimentFromCoolFile(
            file = file,
            resolution = resolution,
            snippets = snippets,
            metadata = metadata,
            topologicalFeatures = topologicalFeatures,
            pairsFile = pairsFile, 
            BPPARAM = BPPARAM
        ))
    }
    if (is_hic(file)) {
        return(.AggrHiCExperimentFromHicFile(
            file = file,
            resolution = resolution,
            snippets = snippets,
            metadata = metadata,
            topologicalFeatures = topologicalFeatures,
            pairsFile = pairsFile, 
            BPPARAM = BPPARAM
        ))
    }
    if (is_hicpro_matrix(file)) {
        if (is.null(bed)) stop("Regions files not provided.")
        if (!is.null(resolution)) stop("Resolution cannot be specified when importing HiC-Pro files.")
        return(.AggrHiCExperimentFromHicproFile(
            file = file,
            bed = bed,
            metadata = metadata,
            snippets = snippets,
            topologicalFeatures = topologicalFeatures,
            pairsFile = pairsFile, 
            BPPARAM = BPPARAM
        ))
    }
    return(NA)
}

#' @rdname AggrHiCExperiment

.AggrHiCExperimentFromCoolFile <- function(
    file, 
    resolution = NULL, 
    snippets, 
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
    BPPARAM
) {
    
    ## -- Check that provided file is valid
    file <- gsub('~', Sys.getenv('HOME'), file)
    check_cool_file(file)
    check_cool_format(file, resolution)

    ## -- Parse multi-query
    if (is(snippets, 'GRanges')) {
        gis <- .coolMulti1DQuery(file, resolution, snippets, BPPARAM = BPPARAM)
        # snippets <- S4Vectors::Pairs(snippets, snippets)
        # gis <- .coolMulti2DQuery(file, resolution, snippets, BPPARAM = BPPARAM)
    }
    else if (is(snippets, 'Pairs')) {
        gis <- .coolMulti2DQuery(file, resolution, snippets, BPPARAM = BPPARAM)
    }
    mcols <- GenomicRanges::mcols(gis)
    GenomicRanges::mcols(gis) <- NULL

    ## -- Create contact object
    x <- methods::new("AggrHiCExperiment", 
        fileName = as.character(file),
        focus = "aggregated snippets", 
        resolutions = lsCoolResolutions(file), 
        resolution = ifelse(is.null(resolution), lsCoolResolutions(file)[1], resolution), 
        interactions = gis, 
        scores = S4Vectors::SimpleList(
            'count' = as.numeric(mcols$count), 
            'balanced' = as.numeric(mcols$balanced), 
            'expected' = as.numeric(mcols$expected), 
            'detrended' = as.numeric(mcols$detrended)
        ),
        topologicalFeatures = c(snippets, topologicalFeatures),
        pairsFile = pairsFile, 
        metadata = metadata
    )
    methods::validObject(x)
    return(x)
} 

#' @rdname AggrHiCExperiment
