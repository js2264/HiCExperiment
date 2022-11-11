#' @importFrom BiocIO import
#' @importFrom stringr str_split
#' @export

setMethod('import', signature = 'CoolFile', function(con, ...) {

    path <- path(con)
    stopifnot(file.exists(path))
    check_cool_format(path, resolution(con))
    HiCExperiment(con, ...)

})

setMethod('import', signature = 'PairsFile', function(con, ...) {

    con <- path(con)
    stopifnot(file.exists(con))
    pairs2gi(con, ...)

})

