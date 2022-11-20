#' @title Checks functions
#' 
#' @name checks
#' @rdname checks
#' 
#' @description 
#' 
#' Useful functions to validate the nature/structure of (m)cool files or 
#' `HiCExperiment` objects.
#'  All these check functions should return a logical.
#' 
#' @param path Path of a (m)cool file
#' @param contacts A `HiCExperiment` object
#' @param resolution Resolution
#' @param pair Pairs object with length of 1
#' @param ... `HiCExperiment` object, arguments passed on by other functions
#' @return Logical
NULL

#' @rdname checks

check_cool_file <- function(path) {
    if (!file.exists(path) | {
        isTRUE(nzchar(Sys.readlink(path), keepNA = TRUE)) & 
        !file.exists(Sys.readlink(path))
    }) {
        stop('File not found. Aborting now')
    }
    if (!{is_cool(path) | is_mcool(path)}) {
        stop('Provided file is not a .cool/.mcool file.\n  Aborting now.')
    }
    TRUE
}

#' @rdname checks

check_resolution <- function(contacts, resolution) {
    available_res <- resolutions(contacts)
    if (!resolution %in% available_res) 
        stop("Resolution not stored in the cool file.\n", paste0('  Available resolutions: ', paste0(res, collapse = ', '), '.'))
    TRUE
}

#' @rdname checks

check_cool_format <- function(path, resolution, ...) {
    if (is_mcool(path)) {
        res <- lsCoolResolutions(path)
        if (is.null(resolution)) {
            stop("File is in .mcool format, a resolution must be provided.\n", paste0('  Available resolutions: ', paste0(res, collapse = ', '), '.'))
        }
        if (!resolution %in% res) {
            stop("Resolution not stored in cool file.\n", paste0('  Available resolutions: ', paste0(res, collapse = ', '), '.'))
        }
    }
    if (is_cool(path)) {
        if (!is.null(resolution)) {
            stop("File is in .cool format, please do not specify any resolution. Aborting now.")
        }
    }
    TRUE
}

#' @rdname checks

is_mcool <- function(path) {
    x <- lsCoolFiles(path)
    all(grepl('^/resolutions', x))
}

#' @rdname checks

is_cool <- function(path) {
    x <- lsCoolFiles(path)
    all(!grepl('^/resolutions', x))
}

#' @rdname checks

is_square <- function(pair) {
    w1 <- GenomicRanges::width(S4Vectors::first(pair))
    w2 <- GenomicRanges::width(S4Vectors::second(pair))
    if (w1 != w2) {
        stop("Provided pair is not square.")
    }
    TRUE
}
