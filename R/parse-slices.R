#' @name multi2Query
#' @title Querying multiple slices of a contact matrix
#' 
#' @description
#' These functions are the workhorse internal functions used to extract 
#' counts from multiple genomic coordinates in a Hi-C contact matrix.
#'
#' @param file path to a Hi-C contact file (can be any format, (m)cool, .hic, or HiC-Pro-derived)
#' @param resolution resolution to use to import matrix over specified targets
#' @param pairs slices to read, provided as a Pairs object
#' @param BPPARAM BiocParallel parameters
#' @param bed associated bed file for HiC-Pro derived contact matrix. 
#' @param maxDistance Maximum distance to use when compiling distance decay
#' @return a GInteractions object with `count`, `balanced`, `detrended` and 
#'   `expected` scores
#'
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel MulticoreParam
#' @keywords internal

.multi2DQuery <- function(
    file, 
    resolution, 
    pairs, 
    maxDistance = NULL, 
    bed = NULL, 
    BPPARAM = BiocParallel::bpparam()  
) {

    message( "Going through preflight checklist..." )
    # - Get si and bins from contact matrix
    if (.is_cool(file)) {
        allBins <- .getCoolAnchors(file, resolution = NULL)
        si <- .cool2seqinfo(file, NULL)
    } else if (.is_mcool(file)) {
        allBins <- .getCoolAnchors(file, resolution = resolution)
        si <- .cool2seqinfo(file, resolution)
    }
    else if (.is_hic(file)) {
        allBins <- .getHicAnchors(file, resolution = resolution)
        si <- .hic2seqinfo(file)
    }
    else if (.is_hicpro_matrix(file) & .is_hicpro_regions(bed)) {
        allBins <- .getHicproAnchors(bed)
        si <- .hicpro2seqinfo(file)
    }
    allBins <- allBins[GenomicRanges::width(allBins) == resolution]
    
    # - Filter and reformat provided pairs
    pairs <- pairs[S4Vectors::queryHits(
        GenomicRanges::findOverlaps(
            S4Vectors::first(pairs), GenomicRanges::reduce(allBins), type = 'within'
        )
    )]
    pairs <- pairs[S4Vectors::queryHits(
        GenomicRanges::findOverlaps(
            S4Vectors::second(pairs), GenomicRanges::reduce(allBins), type = 'within'
        )
    )]
    coordsList1 <- S4Vectors::first(pairs)
    coordsList1$bin_id <- allBins[
        S4Vectors::subjectHits(GenomicRanges::findOverlaps(
            GenomicRanges::resize(coordsList1, fix = 'center', width = 1), allBins))
    ]$bin_id
    coordsList2 <- S4Vectors::second(pairs)
    coordsList2$bin_id <- allBins[
        S4Vectors::subjectHits(GenomicRanges::findOverlaps(
            GenomicRanges::resize(coordsList2, fix = 'center', width = 1), allBins))
    ]$bin_id
    breadth <- GenomicRanges::width(S4Vectors::first(pairs))[1]
    if(is.na(breadth)) stop("No pairs are contained within the contact matrix. Try with a smaller `flankingBins` or a finer `resolution`")
    if (is.null(maxDistance)) {
        maxDistance <- InteractionSet::GInteractions(
            S4Vectors::first(pairs), S4Vectors::second(pairs)
        ) |> InteractionSet::pairdist(type = 'span') |> max()
        if (any(c(is.na(maxDistance), maxDistance < breadth))) maxDistance <- breadth
        maxDistance <- ceiling(maxDistance / resolution) * resolution
    }
    is1D <- all(S4Vectors::first(pairs) == S4Vectors::second(pairs))
    isTrans <- all(GenomicRanges::seqnames(coordsList1) != GenomicRanges::seqnames(coordsList2))

    # - !!! HEAVY LOAD !!! Parse ALL pixels and convert to sparse matrix
    # - Full parsing has to be done since parallelized access to HDF5 is not supported 
    # - Once full parsing is done, parallelization is trivial
    message( "Parsing the entire contact matrice as a sparse matrix..." )
    if (.is_cool(file)) {
        l <- .dumpCool(file, resolution = NULL)
    }
    if (.is_mcool(file)) {
        l <- .dumpCool(file, resolution = resolution)
    }
    else if (.is_hic(file)) {
        l <- .dumpHic(file, resolution = resolution)
    }
    else if (.is_hicpro_matrix(file)) {
        l <- .dumpHicpro(file, bed)
    }

    # - Get detrending model
    message( "Modeling distance decay..." )
    detrendingModel <- distanceDecay(l, maxDistance)
    detrendingModelMat <- .df2symmmat(
        detrendingModel$diag, detrendingModel$score
    )
    detrendingModelMat[lower.tri(detrendingModelMat)] <- NA

    # - For each pair, recover balanced and detrended heatmap
    message( "Filtering for contacts within provided targets..." )
    l_count <- Matrix::sparseMatrix(
        i= l[['pixels']]$bin1_id + 1,
        j= l[['pixels']]$bin2_id + 1,
        x= l[['pixels']]$count
    )
    l_balanced <- Matrix::sparseMatrix(
        i= as.vector(l[['pixels']]$bin1_id + 1),
        j= as.vector(l[['pixels']]$bin2_id + 1),
        x= as.vector(l[['pixels']]$score)
    )
    xx <- BiocParallel::bplapply(BPPARAM = BPPARAM, 
        seq_along(pairs), function(K) {
            pair <- pairs[K]
            dist <- IRanges::start(S4Vectors::second(pair)) - 
                IRanges::start(S4Vectors::first(pair))
            dist <- abs(ceiling(dist / resolution))
            bi_1 <- seq(
                coordsList1[K]$bin_id - {breadth/resolution}/2, 
                coordsList1[K]$bin_id + {breadth/resolution}/2
            )
            bi_2 <- seq(
                coordsList2[K]$bin_id - {breadth/resolution}/2, 
                coordsList2[K]$bin_id + {breadth/resolution}/2
            )
            # bin1_weights <- allBins[bi_1]$weight
            # bin2_weights <- allBins[bi_2]$weight
            # counts <- l_count[bi_1+1, bi_2+1]
            # balanced <- t(apply(
            #     apply(counts, 1, `*`, bin1_weights), 2, `*`, bin2_weights
            # ))
            counts <- l_count[bi_1+1, bi_2+1]
            balanced <- l_balanced[bi_1+1, bi_2+1]
            if (is1D) balanced[lower.tri(balanced)] <- NA
            exp_bi_0 <- seq(
                {{breadth/resolution}/2} - {breadth/resolution}/2, 
                {{breadth/resolution}/2} + {breadth/resolution}/2
            )
            expected <- detrendingModelMat[exp_bi_0+1, exp_bi_0+1+dist]
            if (is1D) expected[lower.tri(expected)] <- NA
            detrended <- log2( 
                {balanced/sum(balanced, na.rm = TRUE)} / 
                {expected/sum(expected, na.rm = TRUE)} 
            )
            if (is1D) detrended[lower.tri(detrended)] <- NA
            detrended[is.infinite(detrended)] <- NA
            list(
                count = counts, 
                balanced = balanced,
                expected = expected, 
                detrended = detrended
            )
        }
    )
    countsMats <- simplify2array(lapply(xx, function(.x) as.matrix(.x[['count']])))
    balancedMats <- simplify2array(lapply(xx, function(.x) as.matrix(.x[['balanced']])))
    expectedMats <- simplify2array(lapply(xx, function(.x) as.matrix(.x[['expected']])))
    detrendedMats <- simplify2array(lapply(xx, function(.x) as.matrix(.x[['detrended']])))
    nNonemptySlices <- sum(unlist(lapply(seq_len(length(pairs)), 
        function(k) {
            sum(countsMats[ , , k], na.rm = TRUE) > 0
        }
    )))
    counts <- apply(countsMats, c(1, 2), sum, na.rm = TRUE) / nNonemptySlices
    balanced <- apply(balancedMats, c(1, 2), sum, na.rm = TRUE) / nNonemptySlices
    expected <- apply(expectedMats, c(1, 2), sum, na.rm = TRUE) / nNonemptySlices
    detrended <- apply(detrendedMats, c(1, 2), sum, na.rm = TRUE) / nNonemptySlices
    if (is1D) counts[lower.tri(counts)] <- NA
    if (is1D) balanced[lower.tri(balanced)] <- NA
    if (is1D) expected[lower.tri(expected)] <- NA
    if (is1D) detrended[lower.tri(detrended)] <- NA

    # - Convert in dense GInteractions
    an <- GenomicRanges::GRanges(
        'aggr.', IRanges::IRanges(
            start = seq(
                -{breadth/resolution}/2,
                +{breadth/resolution}/2, 
                length.out = {breadth/resolution}+1
            ) * resolution + 1 -resolution/2, 
            width = resolution
        ), 
        bin_id = seq(
            -{breadth/resolution}/2, +{breadth/resolution}/2, 
            length.out = {breadth/resolution}+1
        )
    )
    an1 <- rep(an, each = {breadth/resolution}+1)
    an2 <- rep(an, {breadth/resolution}+1)
    gis <- InteractionSet::GInteractions(an1, an2)
    gis$diag <- {IRanges::start(InteractionSet::anchors(gis)[[2]])
        - IRanges::start(InteractionSet::anchors(gis)[[1]])
    } / resolution
    gis$count <- as.vector(counts)
    gis$balanced <- as.vector(balanced)
    gis$expected <- as.vector(expected)
    gis$detrended <- as.vector(detrended)
    gis <- sort(gis)
    S4Vectors::metadata(gis) <- list(slices = S4Vectors::SimpleList(
        count = countsMats,
        balanced = balancedMats,
        expected = expectedMats,
        detrended = detrendedMats 
    ), pairs = pairs)
    return(gis)

}
