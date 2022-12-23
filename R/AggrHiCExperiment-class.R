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

methods::setClass("AggrHiCExperiment", contains = c("HiCExperiment"), 
    slots = c(slices = "SimpleList")
)

#' @rdname AggrHiCExperiment
#' @export
#' @examples
#' library(rtracklayer)
#' devtools::load_all('.')
#' file = '/home/rsg/.cache/R/ExperimentHub/7fa45373d163_7836'
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
#' #snippets <- snippets[strand(snippets) == '+']
#' #snippets <- snippets[resize(snippets, fix = 'center', width = 1) %over% GRanges('chr2:1-12000000')]
#' pairs <- GInteractions(
#'     rep(snippets, each = length(snippets)),
#'     rep(snippets, length(snippets))
#' )
#' pairs <- swapAnchors(pairs)
#' pairs <- pairs[!is.na(pairdist(pairs))]
#' pairs <- pairs[pairdist(pairs) > 30000 & pairdist(pairs) < 500000]

AggrHiCExperiment <- function(
    file, 
    resolution = NULL, 
    snippets, 
    metadata = list(), 
    topologicalFeatures = S4Vectors::SimpleList(), 
    pairsFile = NULL, 
    BPPARAM = BiocParallel::bpparam(),
    ...
) {
    file <- gsub('~', Sys.getenv('HOME'), file)
    stopifnot(file.exists(file))
    params <- list(...)
    bed <- NULL ; if ("bed" %in% names(params)) bed <- params[['bed']]
    if (!is.null(resolution)) resolution <- as.integer(resolution)
    
    if (is_cool(file) | is_mcool(file)) {
        check_cool_file(file)
        check_cool_format(file, resolution)
        return(.aggrHiCExperiment(
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
        check_hic_file(file)
        check_hic_format(file, resolution)
        return(.aggrHiCExperiment(
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
        check_hicpro_files(file, bed)
        if (is.null(bed)) stop("Regions files not provided.")
        return(.aggrHiCExperiment(
            file = file,
            resolution = resolution,
            snippets = snippets,
            metadata = metadata,
            topologicalFeatures = topologicalFeatures,
            pairsFile = pairsFile, 
            bed = bed,
            BPPARAM = BPPARAM
        ))
    }
    return(NA)
}

#' @rdname AggrHiCExperiment

.aggrHiCExperiment <- function(
    file, 
    resolution = NULL, 
    snippets, 
    metadata = list(), 
    topologicalFeatures = S4Vectors::SimpleList(), 
    pairsFile = NULL, 
    BPPARAM
) {
    
    ## -- Parse multi-query
    if (is(snippets, 'GRanges')) {

        # Need to check that snippets are OK (unique width, greater than 0, ...)
        snippets <- S4Vectors::Pairs(snippets, snippets)
        gis <- .multi2DQuery(file, resolution, snippets, BPPARAM = BPPARAM)

    }
    else if (is(snippets, 'GInteractions')) {

        # Need to check that pairs are OK (all width = 1)
        # Need to check that pairs are OK (distance between each pair > width of each anchor)
        gis <- .multi2DQuery(file, resolution, as(snippets, 'Pairs'), BPPARAM = BPPARAM)

    }
    mdata <- S4Vectors::metadata(gis)
    S4Vectors::metadata(gis) <- list()
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
        slices = mdata$slices, 
        topologicalFeatures = S4Vectors::SimpleList(
            c(snippets = snippets, as.list(topologicalFeatures))
        ),
        pairsFile = pairsFile, 
        metadata = metadata
    )
    methods::validObject(x)
    return(x)
} 
