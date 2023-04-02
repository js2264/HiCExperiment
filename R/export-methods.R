
#' @title HiCExperiment export methods
#' 
#' @name export-methods
#' @aliases export
#' @aliases export,HiCExperiment,character,missing-method
#' 
#' @description 
#' 
#' Export methods to save `HiCExperiment` objects into HiC-Pro-like `regions` 
#' and `matrix` files (`regions.bed` and `matrix.mtx`)
#' 
#' @param object HiCExperiment object
#' @param con Prefix of `*_regions.bed` and `*_matrix.mtx` files
#' @usage export(object, con)
#' @return A `HicproFile` object
#' 
#' @importFrom BiocIO export
#' @examples 
#' mcoolPath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' x <- import(mcoolPath, resolution = 16000, focus = 'XVI', format = 'cool')

NULL

#' @exportMethod export 

setMethod('export', signature(object = 'HiCExperiment', con = 'character'), function(object, con) {

    out_regions <-  paste0(con, '_regions.bed')
    out_matrix <-  paste0(con, '_matrix.mtx') 

    # -- Export regions
    b <- bins(object) 
    tab <- data.frame(
        seqnames(b),
        start(b) - 1, 
        end(b), 
        b$bin_id+1
    )
    message("Writing regions to ", out_regions, " file...")
    vroom::vroom_write(
        tab, out_regions, 
        col_names = FALSE, 
        quote = "none", 
        progress = FALSE
    )

    # -- Export matrix
    gi <- interactions(object)
    tab <- data.frame(
        gi$bin_id1+1,
        gi$bin_id2+1,
        gi$count
    )
    message("Writing raw count matrix to ", out_matrix, " file...")
    vroom::vroom_write(
        tab, out_matrix, 
        col_names = FALSE, 
        quote = "none", 
        progress = FALSE
    )

})