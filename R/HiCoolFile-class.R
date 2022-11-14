#' @title `HiCoolFile` S4 class and methods
#'
#' @slot resolution numeric value or NULL 
#' 
#' @importFrom methods setClass
#' @importClassesFrom BiocIO BiocFile
#' @export
#' @rdname HiCoolFile 
#' @include CoolFile-class.R
#' @include CoolFile-method.R
#' @include PairsFile-class.R

setClass('HiCoolFile', contains = c('CoolFile', 'PairsFile'),
    slots = list(
        coolFile = 'CoolFile',
        pairsFile = 'PairsFile',
        log = 'list'
    )
)

#' @rdname HiCoolFile 
#' @export

HiCoolFile <- function(mcool, resolution, pairs, log) {
    new(
        'HiCoolFile', 
        CoolFile(path = mcool, resolution = resolution), 
        pairsFile = PairsFile(pairs), 
        log = log
    )
}
