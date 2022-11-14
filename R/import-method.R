#' @importFrom BiocIO import
#' @importFrom stringr str_split
#' @export

setMethod('import', signature = 'CoolFile', function(con, ...) {

    path <- BiocGenerics::path(con)
    stopifnot(file.exists(path))
    check_cool_format(path, resolution(con))
    HiCExperiment(con, ...)

})

setMethod('import', signature = 'HiCoolFile', function(con, ...) {

    path <- BiocGenerics::path(con)
    stopifnot(file.exists(path))
    check_cool_format(path, resolution(con))
    HiCExperiment(con, ..., pairsFile = BiocGenerics::path(con@pairsFile))

})

setMethod('import', signature = 'PairsFile', function(con, ...) {

    con <- BiocGenerics::path(con)
    stopifnot(file.exists(con))
    pairs2gi(con, ...)

})

