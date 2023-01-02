#' @name multi2Query
#' @title Querying multiple slices of a contact matrix
#' 
#' @description
#' These functions are the workhorse internal functions used to extract 
#' counts from multiple genomic coordinates in a Hi-C contact matrix.
#'
#' @param file file
#' @param resolution resolution
#' @param pairs slices to read
#' @param BPPARAM BiocParallel parameters
#' @param bed associated bed file for HiC-Pro derived contact matrix. 
#' @return a GInteractions object with `count`, `balanced`, `detrended` and 
#'   `expected` scores
#'
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel MulticoreParam

.multi2DQuery <- function(file, resolution, pairs, 
    BPPARAM = BiocParallel::MulticoreParam(workers = 8, progressbar = TRUE), 
    bed = NULL
) {

    message( "Going through preflight checklist..." )
    # - Get si and bins from contact matrix
    if (is_cool(file)) {
        all_bins <- .getCoolAnchors(file, resolution = NULL)
        si <- .cool2seqinfo(file, NULL)
    } else if (is_mcool(file)) {
        all_bins <- .getCoolAnchors(file, resolution = resolution)
        si <- .cool2seqinfo(file, resolution)
    }
    else if (is_hic(file)) {
        all_bins <- .getHicAnchors(file, resolution = resolution)
        si <- .hic2seqinfo(file)
    }
    else if (is_hicpro_matrix(file) & is_hicpro_regions(bed)) {
        all_bins <- .getHicproAnchors(bed)
        si <- .hicpro2seqinfo(file)
    }
    all_bins <- all_bins[GenomicRanges::width(all_bins) == resolution]
    
    # - Filter and reformat provided pairs
    pairs <- pairs[S4Vectors::queryHits(
        GenomicRanges::findOverlaps(
            S4Vectors::first(pairs), GenomicRanges::reduce(all_bins), type = 'within'
        )
    )]
    pairs <- pairs[S4Vectors::queryHits(
        GenomicRanges::findOverlaps(
            S4Vectors::second(pairs), GenomicRanges::reduce(all_bins), type = 'within'
        )
    )]
    coords_list_1 <- S4Vectors::first(pairs)
    coords_list_1$bin_id <- all_bins[
        S4Vectors::subjectHits(GenomicRanges::findOverlaps(
            GenomicRanges::resize(coords_list_1, fix = 'center', width = 1), all_bins))
    ]$bin_id
    coords_list_2 <- S4Vectors::second(pairs)
    coords_list_2$bin_id <- all_bins[
        S4Vectors::subjectHits(GenomicRanges::findOverlaps(
            GenomicRanges::resize(coords_list_2, fix = 'center', width = 1), all_bins))
    ]$bin_id
    breadth <- GenomicRanges::width(S4Vectors::first(pairs))[1]
    if(is.na(breadth)) stop("No pairs are contained within the contact matrix. Try with a smaller `flanking_bins` or a finer `resolution`")
    threshold <- InteractionSet::GInteractions(
        S4Vectors::first(pairs), S4Vectors::second(pairs)
    ) |> InteractionSet::pairdist(type = 'span') |> max()
    if (threshold < breadth) threshold <- breadth
    threshold <- ceiling(threshold / resolution) * resolution
    is1D <- all(S4Vectors::first(pairs) == S4Vectors::second(pairs))

    # - !!! HEAVY LOAD !!! Parse ALL pixels and convert to sparse matrix
    # - Full parsing has to be done since parallelized access to HDF5 is not supported 
    # - Once full parsing is done, parallelization is trivial
    message( "Parsing the entire contact matrice as a sparse matrix..." )
    if (is_cool(file)) {
        l <- .dumpCool(file, resolution = NULL)
    }
    if (is_mcool(file)) {
        l <- .dumpCool(file, resolution = resolution)
    }
    else if (is_hic(file)) {
        l <- .dumpHic(file, resolution = resolution)
    }
    else if (is_hicpro_matrix(file)) {
        l <- .dumpHicpro(file, bed)
    }

    # - Get detrending model
    message( "Modeling distance decay..." )
    detrending_model <- distance_decay(l, threshold)
    detrending_model_mat <- .df_to_symmmat(
        detrending_model$diag, detrending_model$score
    )
    detrending_model_mat[lower.tri(detrending_model_mat)] <- NA

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
                coords_list_1[K]$bin_id - {breadth/resolution}/2, 
                coords_list_1[K]$bin_id + {breadth/resolution}/2
            )
            bi_2 <- seq(
                coords_list_2[K]$bin_id - {breadth/resolution}/2, 
                coords_list_2[K]$bin_id + {breadth/resolution}/2
            )
            # bin1_weights <- all_bins[bi_1]$weight
            # bin2_weights <- all_bins[bi_2]$weight
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
            expected <- detrending_model_mat[exp_bi_0+1, exp_bi_0+1+dist]
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
    counts_mats <- simplify2array(lapply(xx, function(.x) as.matrix(.x[['count']])))
    balanced_mats <- simplify2array(lapply(xx, function(.x) as.matrix(.x[['balanced']])))
    expected_mats <- simplify2array(lapply(xx, function(.x) as.matrix(.x[['expected']])))
    detrended_mats <- simplify2array(lapply(xx, function(.x) as.matrix(.x[['detrended']])))
    n_nonempty_slices <- sum(unlist(lapply(seq_len(length(pairs)), 
        function(k) {
            sum(counts_mats[ , , k], na.rm = TRUE) > 0
        }
    )))
    counts <- apply(counts_mats, c(1, 2), sum, na.rm = TRUE) / n_nonempty_slices
    balanced <- apply(balanced_mats, c(1, 2), sum, na.rm = TRUE) / n_nonempty_slices
    expected <- apply(expected_mats, c(1, 2), sum, na.rm = TRUE) / n_nonempty_slices
    detrended <- apply(detrended_mats, c(1, 2), sum, na.rm = TRUE) / n_nonempty_slices
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
    an1 <- rep(an, {breadth/resolution}+1)
    an2 <- rep(an, each = {breadth/resolution}+1)
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
        count = counts_mats,
        balanced = balanced_mats,
        expected = expected_mats,
        detrended = detrended_mats 
    ))
    return(gis)

}
