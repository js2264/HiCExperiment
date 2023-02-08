#' @rdname import-methods
#' @export 

setMethod("availableResolutions", signature(x = "ANY"), function(x, ...) {
    .availableResolutions(x, ...)
})

#' @rdname import-methods
#' @export 
#' @include CoolFile-class.R
#' @include import-methods.R

setMethod("availableResolutions", signature(x = "CoolFile"), function(x) {
    .availableResolutions(path(x))
})

#' @rdname import-methods
#' @export 
#' @include HicFile-class.R

setMethod("availableResolutions", signature(x = "HicFile"), function(x) {
    .availableResolutions(path(x))
})

#' @rdname import-methods
#' @export 
#' @include HicproFile-class.R

setMethod("availableResolutions", signature(x = "HicproFile"), function(x) {
    .availableResolutions(path(x), x@bed)
})


#' @rdname import-methods
#' @export 

setMethod("availableChromosomes", signature(x = "ANY"), function(x, ...) {
    .availableChromosomes(x, ...)
})

#' @rdname import-methods
#' @export 
#' @include CoolFile-class.R

setMethod("availableChromosomes", signature(x = "CoolFile"), function(x) {
    .availableChromosomes(path(x))
})

#' @rdname import-methods
#' @export 
#' @include HicFile-class.R

setMethod("availableChromosomes", signature(x = "HicFile"), function(x) {
    .availableChromosomes(path(x))
})

#' @rdname import-methods
#' @export 
#' @include HicproFile-class.R

setMethod("availableChromosomes", signature(x = "HicproFile"), function(x) {
    .availableChromosomes(path(x), x@bed)
})


.availableResolutions <- function(file, bed = NULL) {
    file <- gsub('~', Sys.getenv('HOME'), file)
    stopifnot(file.exists(file))
    if (is_cool(file) | is_mcool(file)) {
        return(lsCoolResolutions(file, verbose = TRUE))
    }
    if (is_hic(file)) {
        return(lsHicResolutions(file, verbose = TRUE))
    }
    if (is_hicpro_matrix(file)) {
        if (is.null(bed)) stop("Regions files not provided.")
        return({
            bed1 <- vroom::vroom(
                bed, 
                col_names = FALSE, 
                progress = FALSE, 
                show_col_types = FALSE, 
                n_max = 10
            )
            max(unique(bed1[[3]][1] - bed1[[2]][1]))
        })
    }
}

.availableChromosomes <- function(file, bed = NULL) {
    file <- gsub('~', Sys.getenv('HOME'), file)
    stopifnot(file.exists(file))
    if (is_cool(file)) {
        return(.cool2seqinfo(file))
    }
    if (is_mcool(file)) {
        return(.cool2seqinfo(file, resolution = lsCoolResolutions(file)[1]))
    }
    if (is_hic(file)) {
        return(.hic2seqinfo(file))
    }
    if (is_hicpro_matrix(file)) {
        if (is.null(bed)) stop("Regions files not provided.")
        return(.hicpro2seqinfo(bed))
    }
}

