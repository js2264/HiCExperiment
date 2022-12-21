.coolMulti1DQuery <- function(file, resolution, snippets, BPPARAM = BiocParallel::bpparam()) {

    message( "Going through preflight checklist..." )
    # - Get bins from cool
    all_bins <- .getCoolAnchors(file, resolution = resolution)
    all_bins <- all_bins[width(all_bins) == resolution]
    
    # - Filter and reformat provided snippets
    si <- .cool2seqinfo(file, resolution)
    coords_list <- subsetByOverlaps(snippets, reduce(all_bins), type = 'within')
    GenomeInfoDb::seqlevels(coords_list) <- GenomeInfoDb::seqlevels(si)
    coords_list <- sort(coords_list)
    coords_list$bin_id <- all_bins[
        S4Vectors::subjectHits(GenomicRanges::findOverlaps(
            resize(coords_list, fix = 'center', width = 1), all_bins))
    ]$bin_id
    coords_list <- sort(coords_list)
    breadth <- width(coords_list)[1]
    
    # - !!! HEAVY LOAD !!! Parse ALL pixels and filter those corresponding to snippets
    message( "Parsing the entire contact matrice as a sparse matrix..." )
    l <- .dumpCool(file, resolution = resolution)
    l <- Matrix::sparseMatrix(
        i= l[['pixels']]$bin1_id + 1,
        j= l[['pixels']]$bin2_id + 1,
        x= l[['pixels']]$count
    )

    # - For each snippet, recover corresponding pixels
    message( "Filtering for contacts within provided snippets..." )
    xx <- BiocParallel::bplapply(BPPARAM = BPPARAM, 
        seq_along(coords_list), function(K) {
            bi_1 <- seq(
                coords_list[K]$bin_id - {breadth/resolution}/2, 
                coords_list[K]$bin_id + {breadth/resolution}/2
            )
            bi_2 <- bi_1
            nbins <- breadth/resolution+1
            bin1_weights <- all_bins[bi_1+1]$weight
            bin2_weights <- all_bins[bi_2+1]$weight
            counts <- l[bi_1+1, bi_2+1]
            counts[lower.tri(counts)] <- NA
            balanced <- t(apply(
                apply(counts, 1, `*`, bin1_weights), 2, `*`, bin2_weights
            ))
            list(
                count = counts,
                balanced = balanced
            )
        }
    )
    counts <- apply(simplify2array(
        purrr::map(xx, ~as.matrix(.x[['count']]))), 
        c(1, 2), sum, na.rm = TRUE) / length(coords_list
    )
    balanced <- apply(simplify2array(
        purrr::map(xx, ~as.matrix(.x[['balanced']]))), 
        c(1, 2), sum, na.rm = TRUE) / length(coords_list
    )
    
    # - Convert in dense GInteractions
    an <- GenomicRanges::GRanges(
        'aggr.', IRanges::IRanges(
            start = seq(-{breadth/resolution}/2, +{breadth/resolution}/2, length.out = {breadth/resolution}+1) * resolution + 1, 
            width = resolution
        ), 
        bin_id = seq(-{breadth/resolution}/2, +{breadth/resolution}/2, length.out = {breadth/resolution}+1)
    )
    an1 <- rep(an, {breadth/resolution}+1)
    an2 <- rep(an, each = {breadth/resolution}+1)
    gis <- InteractionSet::GInteractions(an1, an2)
    gis$diag <- {start(anchors(gis)[[2]]) - start(anchors(gis)[[1]])} / resolution
    gis$count <- as.vector(counts)
    gis$balanced <- as.vector(balanced)
    gis <- sort(gis[gis$diag >= 0])

    # - Detrend GInteractions
    expected <- tibble::as_tibble(gis) |> 
        dplyr::group_by(diag) |>
        dplyr::summarize(average_interaction_per_diag = mean(balanced, na.rm = TRUE)) |> 
        dplyr::mutate(average_interaction_per_diag = average_interaction_per_diag / 2)
    gis$expected <- tibble::as_tibble(gis) |> 
        dplyr::left_join(expected, by = 'diag') |> 
        dplyr::pull(average_interaction_per_diag)
    gis$detrended <- log2( {gis$balanced/sum(gis$balanced)} / {gis$expected/sum(gis$expected)} )

    return(gis)

}




