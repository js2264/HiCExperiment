#' @name HiCExperiment_example
#'
#' @title Manage cache / download files for HiCExperiment package
#'
#' @description Managing data downloads is important to save disk space and
#' re-downloading data files. This can be done effortlessly via the integrated
#' `BiocFileCache` system.
#'
#' @section HiCExperimentCache:
#' Get the directory location of the cache. It will prompt the user to create
#' a cache if not already created. A specific directory can be used via
#' \code{setCache}.
#'
#' @section setCache:
#' Specify the directory location of the data cache. By default, it will
#' go into the user's home and package name directory as given by
#' \link[tools:userdir]{R_user_dir} (default: varies by system e.g., for Linux:
#' '$HOME/.cache/R/HiCExperiment').
#'
#' @param directory character(1) The file location where the cache is located.
#' Once set, future downloads will go to this folder. See `setCache` section
#' for details.
#'
#' @param ask logical(1) (default TRUE when `interactive()`) Confirm the file
#' location of the cache directory
#' 
#' @examples
#' getOption("HiCExperimentCache")
#' HiCExperimentCache()
#'
#' @return The directory / option of the cache location
#'
#' @importFrom GenomicRanges resize
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges 
#' @importFrom IRanges overlapsAny
#' @importFrom InteractionSet GInteractions
#' @importFrom InteractionSet swapAnchors
#' @importFrom InteractionSet pairdist
#' @importFrom S4Vectors SimpleList
#' @importFrom tools R_user_dir
#' @export

HiCExperimentCache <- function(
    directory= tools::R_user_dir("HiCExperiment", "cache")
) {
    cache <- getOption("HiCExperimentCache", .setHiCExperimentCache(directory))
    BiocFileCache::BiocFileCache(cache)
}

#' @rdname HiCExperiment_example

.setHiCExperimentCache <- function(directory, ask = interactive()) {
    stopifnot(
        is.character(directory), length(directory) == 1L, !is.na(directory)
    )
    if (!dir.exists(directory)) {
        if (ask) {
            qtxt <- sprintf(
                "Create HiCExperiment cache at \n    %s? [y/n]: ",
                directory
            )
            if (interactive()) {
                repeat {
                    cat(qtxt)
                    answer <- readLines(n = 1)
                    if (answer %in% c('y', 'Y', 'n', 'N'))
                        break
                }
                tolower(answer)
            } else {
                "n"
            }
            if ("n" == answer)
                stop("'HiCExperiment' directory not created. Use 'setCache'")
        }
        dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    }
    options("HiCExperiment" = directory)
    invisible(directory)
}

#' @rdname HiCExperiment_example
#' @export

HiCExperiment_example <- function() {
    if (!requireNamespace("rtracklayer")) 
        stop("rtracklayer package must be installed.")
    if (!requireNamespace("BiocFileCache")) 
        stop("BiocFileCache package must be installed.")

    bfc <- HiCExperimentCache()

    url_map <- 'https://osf.io/3h9js/download'
    url_ctcf_peaks <- 'https://osf.io/c9pwe/download'
    url_ctcf_track <- 'https://osf.io/w92u3/download'

    # - Fetch contact map and CTCF files
    rid_map <- bfcquery(bfc, query = 'microC.mcool')$rid
    if (!length(rid_map)) {
        message( "Fetching Hi-C contact map from OSF portal" )
        bfcentry <- bfcadd( 
            bfc, 
            rname = 'microC.mcool', 
            fpath = url_map 
        )
        rid_map <- names(bfcentry)
    } 
    rid_ctcf_peaks <- bfcquery(bfc, query = 'CTCF_peaks.narrowPeak')$rid
    if (!length(rid_ctcf_peaks)) {
        message( "Fetching CTCF peaks from OSF portal" )
        bfcentry <- bfcadd( 
            bfc, 
            rname = 'CTCF_peaks.narrowPeak', 
            fpath = url_ctcf_peaks 
        )
        rid_ctcf_peaks <- names(bfcentry)
    }
    rid_ctcf_track <- bfcquery(bfc, query = 'CTCF_track.bw')$rid
    if (!length(rid_ctcf_track)) {
        message( "Fetching CTCF peaks from OSF portal" )
        bfcentry <- bfcadd( 
            bfc, 
            rname = 'CTCF_track.bw', 
            fpath = url_ctcf_track 
        )
        rid_ctcf_track <- names(bfcentry)
    }

    # - Import all files in memory
    message( "Importing CTCF peaks in memory..." )
    res <- 10000
    CTCF_peaks <- rtracklayer::import(
        bfcrpath(bfc, rids = rid_ctcf_peaks), 
        format = 'bed'
    )
    message( "Importing CTCF track in memory..." )
    CTCF_track <- rtracklayer::import(
        bfcrpath(bfc, rids = rid_ctcf_track), 
        format = 'BigWig', 
        selection = rtracklayer::BigWigSelection(
            GenomicRanges::resize(CTCF_peaks, fix = 'center', width = 500)
        ), 
        as = 'NumericList'
    )
    CTCF_peaks$center <- GenomicRanges::start(CTCF_peaks) + {
        GenomicRanges::end(CTCF_peaks) - GenomicRanges::start(CTCF_peaks)
    }/2
    CTCF_peaks$ChIP_score <- BiocGenerics::mean(CTCF_track)
    CTCF_peaks <- GenomicRanges::resize(CTCF_peaks, 
        fix = 'center', width = 200000
    )
    CTCF_peaks <- CTCF_peaks[
        CTCF_peaks$score >= stats::quantile(CTCF_peaks$score , 0.75) & 
        CTCF_peaks$ChIP_score >= stats::quantile(CTCF_peaks$ChIP_score , 0.75)
    ]
    CTCF_peaks <- sort(CTCF_peaks)
    CTCF_peaks <- CTCF_peaks[
        IRanges::overlapsAny(
            GenomicRanges::resize(CTCF_peaks, fix = 'center', width = 1),
            GenomicRanges::GRanges(
                c('chr2', 'chr17'), IRanges::IRanges(
                    start = c(1L, 1L), width = c(242193529L, 83257441L)
                )
            ), 
            type = 'within'
        )
    ]
    CTCF_pairs <- InteractionSet::GInteractions(
        rep(CTCF_peaks, each = length(CTCF_peaks)),
        rep(CTCF_peaks, length(CTCF_peaks))
    )
    CTCF_pairs <- InteractionSet::swapAnchors(CTCF_pairs)
    CTCF_pairs <- CTCF_pairs[!duplicated(CTCF_pairs)]
    CTCF_pairs <- CTCF_pairs[!is.na(InteractionSet::pairdist(CTCF_pairs))]
    CTCF_pairs <- CTCF_pairs[
        InteractionSet::pairdist(CTCF_pairs) > 230000 & 
        InteractionSet::pairdist(CTCF_pairs) < 500000
    ]

    message( "Importing contacts in memory..." )
    x <- HiCExperiment(
        bfcrpath(bfc, rids = rid_map), 
        resolution = res, 
        metadata = list(CTCF_track = bfcrpath(bfc, rids = rid_ctcf_track)), 
        topologicalFeatures = S4Vectors::SimpleList(
            CTCF_peaks = CTCF_peaks, 
            CTCF_pairs = CTCF_pairs
        )
    )

    return(x)
}
