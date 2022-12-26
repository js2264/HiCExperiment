.multi2DQuery <- function(file, resolution, pairs, BPPARAM = BiocParallel::MulticoreParam(workers = 8, progressbar = TRUE), ...) {

    message( "Going through preflight checklist..." )
    # - Get si and bins from contact matrix
    if (is_cool(file) | is_mcool(file)) {
        all_bins <- .getCoolAnchors(file, resolution = resolution)
        si <- .cool2seqinfo(file, resolution)
    }
    else if (is_hic(file)) {
        all_bins <- .getHicAnchors(file, resolution = resolution)
        si <- .hic2seqinfo(file)
    }
    else if (is_hicpro_matrix(file) & is_hicpro_regions(...)) {
        all_bins <- .getHicproAnchors(...)
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
    threshold <- InteractionSet::GInteractions(
        S4Vectors::first(pairs), S4Vectors::second(pairs)
    ) |> InteractionSet::pairdist(type = 'span') |> max()
    if (threshold < breadth) threshold <- breadth
    threshold <- ceiling(threshold / resolution) * resolution

    # - !!! HEAVY LOAD !!! Parse ALL pixels and convert to sparse matrix
    # - Full parsing has to be done since parallelized access to HDF5 is not supported 
    # - Once full parsing is done, parallelization is trivial
    message( "Parsing the entire contact matrice as a sparse matrix..." )
    if (is_cool(file) | is_mcool(file)) {
        l <- .dumpCool(file, resolution = resolution)
    }
    else if (is_hic(file)) {
        l <- .dumpHic(file, resolution = resolution)
    }
    else if (is_hic_pro(file)) {
        l <- .dumpHicpro(file, ...)
    }

    # - Get detrending model
    message( "Modeling distance decay..." )
    detrending_model <- distance_decay(l, threshold)
    detrending_model_mat <- .df_to_symmmat(
        detrending_model$diag, detrending_model$score
    )
    detrending_model_mat[lower.tri(detrending_model_mat)] <- NA

    # - For each pair, recover balanced pixels
    message( "Filtering for contacts within provided pairs..." )
    l <- Matrix::sparseMatrix(
        i= l[['pixels']]$bin1_id + 1,
        j= l[['pixels']]$bin2_id + 1,
        x= l[['pixels']]$count
    )
    xx <- BiocParallel::bplapply(BPPARAM = BPPARAM, 
        seq_along(pairs), function(K) {
            pair <- pairs[K]
            dist <- IRanges::start(S4Vectors::second(pair)) - IRanges::start(S4Vectors::first(pair))
            dist <- abs(ceiling(dist / resolution))
            bi_1 <- seq(
                coords_list_1[K]$bin_id - {breadth/resolution}/2, 
                coords_list_1[K]$bin_id + {breadth/resolution}/2
            )
            bi_2 <- seq(
                coords_list_2[K]$bin_id - {breadth/resolution}/2, 
                coords_list_2[K]$bin_id + {breadth/resolution}/2
            )
            bin1_weights <- all_bins[bi_1]$weight
            bin2_weights <- all_bins[bi_2]$weight
            counts <- l[bi_1+1, bi_2+1]
            balanced <- t(apply(
                apply(counts, 1, `*`, bin1_weights), 2, `*`, bin2_weights
            ))
            balanced[lower.tri(balanced)] <- NA
            exp_bi_0 <- seq(
                {{breadth/resolution}/2} - {breadth/resolution}/2, 
                {{breadth/resolution}/2} + {breadth/resolution}/2
            )
            expected <- detrending_model_mat[exp_bi_0+1, exp_bi_0+1+dist]
            expected[lower.tri(expected)] <- NA
            detrended <- log2( 
                {balanced/sum(balanced, na.rm = TRUE)} / 
                {expected/sum(expected, na.rm = TRUE)} 
            )
            detrended[lower.tri(detrended)] <- NA
            detrended[is.infinite(detrended)] <- NA
            list(
                count = counts, 
                balanced = balanced,
                expected = expected, 
                detrended = detrended
            )
        }
    )
    counts <- apply(simplify2array(
        purrr::map(xx, ~as.matrix(.x[['count']]))
    ), c(1, 2), sum, na.rm = TRUE) / length(pairs)
    balanced <- apply(simplify2array(
        purrr::map(xx, ~as.matrix(.x[['balanced']]))
    ), c(1, 2), sum, na.rm = TRUE) / length(pairs)
    expected <- apply(simplify2array(
        purrr::map(xx, ~as.matrix(.x[['expected']]))
    ), c(1, 2), sum, na.rm = TRUE) / length(pairs)
    detrended <- apply(simplify2array(
        purrr::map(xx, ~as.matrix(.x[['detrended']]))
    ), c(1, 2), sum, na.rm = TRUE) / length(pairs)
    counts[lower.tri(counts)] <- NA
    balanced[lower.tri(balanced)] <- NA
    expected[lower.tri(expected)] <- NA
    detrended[lower.tri(detrended)] <- NA

    # - Convert in dense GInteractions
    an <- GenomicRanges::GRanges(
        'aggr.', IRanges::IRanges(
            start = seq(
                -{breadth/resolution}/2,
                +{breadth/resolution}/2, 
                length.out = {breadth/resolution}+1
            ) * resolution + 1, 
            width = resolution
        ), 
        bin_id = seq(-{breadth/resolution}/2, +{breadth/resolution}/2, length.out = {breadth/resolution}+1)
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
        count = simplify2array(purrr::map(xx, ~as.matrix(.x[['count']]))),
        balanced = simplify2array(purrr::map(xx, ~as.matrix(.x[['balanced']]))),
        expected = simplify2array(purrr::map(xx, ~as.matrix(.x[['expected']]))),
        detrended = simplify2array(purrr::map(xx, ~as.matrix(.x[['detrended']])))
    ))
    return(gis)

}
