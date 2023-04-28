#' @title HiCExperiment export methods
#' 
#' @name export-methods
#' @rdname export-methods
#' @aliases export
#' @aliases export,HiCExperiment,missing,character-method
#' @description 
#' 
#' Export methods to save a HiCExperiment object into a set of HiC-Pro-style 
#' files (matrix & regions files)
#' @param object A HiCExperiment object 
#' @param prefix Prefix used when generating output file(s).
#' @param format File format. Available: `cool` and `HiC-Pro`.
#' @param ... Extra arguments to use when exporting to `cool`. 
#' Can be `metadata <string>` or `chunksize <integer>`.
#' @return Path to saved files
#' @importFrom BiocIO export
#' @examples
#' ################################################################
#' ## ----------- Importing .(m)cool contact matrices ---------- ##
#' ################################################################
#' 
#' mcoolPath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' hic <- import(mcoolPath, format = 'mcool', resolution = 16000)
#' export(hic["II"], prefix = 'subset_chrII', format = 'cool')
#' export(hic["II"], prefix = 'subset_chrII', format = 'HiC-Pro')
NULL

#' @exportMethod export 
#' @rdname export-methods

setMethod("export", signature(object = "HiCExperiment", con = "missing", format = "character"),
    function(object, prefix, format, ...) {

        params <- list(...)

        if (format == 'HiC-Pro') {

            out_regions <- paste0(prefix, "_regions.bed")
            out_matrix <- paste0(prefix, "_matrix.mtx")
            b <- bins(object)
            gi <- interactions(object)
            .f <- .writeHicpro(b, gi, out_regions, out_matrix)

            ## -- Return HicproFile object
            res <- HicproFile(path = out_matrix, bed = out_regions)
            return(res)

        }

        else if (format == 'cool') {

            out_cool <- paste0(prefix, ".cool")
            re <- regions(object)
            b <- as.data.frame(bins(object))
            gi <- interactions(object)

            ## -- Estimate chunksize for pixels
            if ('chunksize' %in% names(params)) {
                chunksize <- params[['chunksize']]
            } else {
                chunksize <- ceiling(sqrt(length(gi)))
                chunksize <- max(1000, chunksize)
            }

            ## -- Get metadata
            if ('metadata' %in% names(params)) {
                metadata <- params[['metadata']]
            } else {
                metadata <- "{}"
            }

            ## -- Create info 
            info <- list(
                "bin-size"= resolution(object),
                "bin-type"= "fixed",
                "creation-date"= gsub(' ', 'T', as.character(Sys.time())),
                "format"= "HDF5::Cooler",
                "format-url"= "https://github.com/open2c/cooler",
                "format-version"= "3",
                "generated-by"= paste0('HiCExperiment-', packageVersion('HiCExperiment')),
                "genome-assembly"= "unknown",
                "metadata" = metadata,
                "nbins"= nrow(b),
                "nchroms"= length(GenomeInfoDb::seqlevels(re)),
                "nnz"= length(interactions(object)),
                "storage-mode"= "symmetric-upper",
                "sum"= sum(scores(object, 'count'))
            )

            ## -- Create Chroms table
            chroms <- data.frame(
                seqnames = GenomeInfoDb::seqlevels(re), 
                seqlengths = GenomeInfoDb::seqlengths(re)
            )
            colnames(chroms) <- c("name", "length")

            ## -- Create Bins table
            bins <- data.frame(
                chrom = as.character(b$seqnames), 
                start = as.integer(b$start - 1),
                end = as.integer(b$end)
            )
            if ('weight' %in% colnames(b)) {
                bins$weight <- b$weight
            }
            else {
                bins$weight <- 1
            }
            bins <- bins[, c("chrom", "end", "start", "weight")]

            ## -- Create Pixels table
            pixels <- data.frame(gi$bin_id1, gi$bin_id2, gi$count)
            colnames(pixels) <- c("bin1_id", "bin2_id", "count")

            ## -- Aggregate all data 
            .f <- .writeCool(info, bins, chroms, pixels, chunksize, out_cool)

            ## -- Return CoolFile object
            res <- CoolFile(.f)
            return(res)

        }

    }
)

.writeHicpro <- function(bins, interactions, out_regions, out_matrix) {
    tab <- data.frame(
        seqnames(bins), start(bins) - 1, end(bins),
        bins$bin_id + 1
    )
    message("Writing regions to ", out_regions, " file...")
    vroom::vroom_write(
        tab, out_regions, 
        col_names = FALSE,
        quote = "none", 
        progress = FALSE
    )
    tab <- data.frame(interactions$bin_id1, interactions$bin_id2, interactions$count)
    message("Writing raw count matrix to ", out_matrix, " file...")
    vroom::vroom_write(
        tab, 
        out_matrix, 
        col_names = FALSE,
        quote = "none", 
        progress = FALSE
    )
}

#' @import rhdf5 

.writeCool <- function(info, bins, chroms, pixels, chunksize, output) {

    unlink(output)
    rhdf5::h5closeAll()
    fid <- rhdf5::H5Fcreate(output)

    ## -- Saving Attributes
    gid <- H5Gopen(fid, "/")
    attrs <- names(info)
    for (attr in attrs) {
        if (is.numeric(info[[attr]])) {
            aid <- rhdf5::H5Acreate(
                h5obj = fid, 
                name = attr, 
                dtype_id = rhdf5::H5Tcopy("H5T_STD_I64LE"), 
                h5space = rhdf5::H5Screate(type = rhdf5::h5default("H5S"), native = TRUE)
            )
            rhdf5::H5Awrite(aid, info[[attr]])
            rhdf5::H5Aclose(aid)
        }
        else {
            rhdf5::h5writeAttribute(
                info[[attr]], 
                gid, 
                attr, 
                encoding="UTF8", 
                variableLengthString=TRUE, 
                asScalar=TRUE
            )
        }
    }
    rhdf5::H5Gclose(gid)

    ## -- Saving bins
    message("Writing /bins table to ", output, " file...")
    gid <- rhdf5::H5Gcreate(fid, 'bins')
    #  chrom
    tid <- rhdf5::H5Tenum_create(dtype_id = "H5T_NATIVE_UCHAR")
    for (i in seq_along(chroms$name)) {
        rhdf5::H5Tenum_insert(
            dtype_id = tid, 
            name = chroms$name[i], 
            value = as.integer(i-1)
        )
    }
    did <- rhdf5::H5Dcreate(
        h5loc = fid, 
        name = "bins/chrom", 
        dtype_id = tid, 
        h5space = rhdf5::H5Screate_simple(dims = length(bins$chrom), native = TRUE)
    )
    rhdf5::H5Dwrite(
        did, 
        as.raw(as.integer(factor(bins$chrom, levels = chroms$name))-1), 
        h5type = tid
    )
    rhdf5::H5Dclose(did)
    #  end
    did <- H5Dcreate(
        h5loc = fid, 
        name = "bins/end", 
        dtype_id = "H5T_NATIVE_INT", 
        h5space = rhdf5::H5Screate_simple(length(bins$end))
    )
    rhdf5::H5Dwrite(
        did, 
        as.integer(bins$end)
    )
    rhdf5::H5Dclose(did)
    #  start
    did <- rhdf5::H5Dcreate(
        h5loc = fid, 
        name = "bins/start", 
        dtype_id = "H5T_NATIVE_INT", 
        h5space = rhdf5::H5Screate_simple(length(bins$start))
    )
    rhdf5::H5Dwrite(
        did, 
        as.integer(bins$start)
    )
    rhdf5::H5Dclose(did)
    #  weight
    #   attributes 
    # attrs <- c("dim", "cis_only", "converged", "ignore_diags", "mad_max", "min_count", "min_nnz", "scale", "tol", "var")
    # for (attr in attrs) {
    #     if (is.numeric(info[[attr]])) {
    #         aid <- rhdf5::H5Acreate(
    #             h5obj = fid, 
    #             name = attr, 
    #             dtype_id = rhdf5::H5Tcopy("H5T_STD_I64LE"), 
    #             h5space = rhdf5::H5Screate(type = h5default("H5S"), native = TRUE)
    #         )
    #         rhdf5::H5Awrite(aid, info[[attr]])
    #         rhdf5::H5Aclose(aid)
    #     }
    #     else {
    #         rhdf5::h5writeAttribute(
    #             info[[attr]], 
    #             gid, 
    #             attr, 
    #             encoding="UTF8", 
    #             variableLengthString=TRUE, 
    #             asScalar=TRUE
    #         )
    #     }
    # }
    #   data
    did <- rhdf5::H5Dcreate(
        h5loc = fid, 
        name = "bins/weight", 
        dtype_id = "H5T_IEEE_F64LE", 
        h5space = rhdf5::H5Screate_simple(length(bins$weight))
    )
    rhdf5::H5Dwrite(
        did, 
        bins$weight
    )
    rhdf5::H5Dclose(did)
    rhdf5::H5Gclose(gid)

    ## -- Saving chroms
    message("Writing /chroms table to ", output, " file...")
    gid <- rhdf5::H5Gcreate(fid, 'chroms')
    #  length
    did <- rhdf5::H5Dcreate(
        h5loc = fid, 
        name = "chroms/length", 
        dtype_id = "H5T_STD_I32LE", 
        h5space = H5Screate_simple(length(chroms$length))
    )
    rhdf5::H5Dwrite(
        did, 
        as.integer(chroms$length)
    )
    rhdf5::H5Dclose(did)
    #  name 
    tid <- rhdf5::H5Tcopy(dtype_id = "H5T_C_S1")
    rhdf5::H5Tset_strpad(tid, strpad = "NULLPAD")
    rhdf5::H5Tset_size(tid, size = max(nchar(chroms$name)))
    rhdf5::H5Tset_cset(tid, cset = 'ASCII')
    did <- rhdf5::H5Dcreate(
        h5loc = fid, 
        name = "chroms/name", 
        dtype_id = tid, 
        h5space = rhdf5::H5Screate_simple(length(chroms$name))
    )
    rhdf5::H5Dwrite(
        did, 
        chroms$name
    )
    rhdf5::H5Dclose(did)
    rhdf5::H5Gclose(gid)

    ## -- Saving pixels
    message("Writing /pixels table to ", output, " file...")
    gid <- rhdf5::H5Gcreate(fid, 'pixels')
    #  bin1_id
    sid <- rhdf5::H5Screate_simple(
        dims = nrow(pixels),
        maxdims = nrow(bins)*ceiling(nrow(bins)/2), 
        native = FALSE
    )
    tid <- "H5T_STD_I64LE"
    rhdf5::h5createDataset(
        file = gid, 
        dataset = 'bin1_id', 
        dims = nrow(pixels), 
        maxdims = nrow(bins)*ceiling(nrow(bins)/2), 
        H5type = tid, 
        chunk = c(chunksize), 
        filter = "DEFLATE", 
        shuffle = TRUE
    )
    rhdf5::h5write(pixels$bin1_id, gid, "bin1_id")
    rhdf5::h5createDataset(
        file = gid, 
        dataset = 'bin2_id', 
        dims = nrow(pixels), 
        maxdims = nrow(bins)*ceiling(nrow(bins)/2), 
        H5type = tid, 
        chunk = c(chunksize), 
        filter = "DEFLATE", 
        shuffle = TRUE
    )
    rhdf5::h5write(pixels$bin2_id, gid, "bin2_id")
    rhdf5::h5createDataset(
        file = gid, 
        dataset = 'count', 
        dims = nrow(pixels), 
        maxdims = nrow(bins)*ceiling(nrow(bins)/2), 
        H5type = "H5T_STD_I32LE", 
        chunk = c(chunksize), 
        filter = "DEFLATE", 
        shuffle = TRUE
    )
    rhdf5::h5write(pixels$count, gid, "count")
    rhdf5::H5Sclose(sid)
    rhdf5::H5Gclose(gid)

    ## -- Saving indexes
    message("Writing /indexes table to ", output, " file...")
    gid <- rhdf5::H5Gcreate(fid, 'indexes')
    indices <- list(
        chrom_offset = c(
            0, cumsum(S4Vectors::runLength(S4Vectors::Rle(bins$chrom)))
        ),
        bin1_offset = c(
            0, factor(pixels$bin1_id+1, levels = seq_len(nrow(bins))) |> 
                table() |> as.vector() |> cumsum()
        )
    )
    #  bin1_offset 
    did <- rhdf5::H5Dcreate(
        h5loc = fid, 
        name = "indexes/bin1_offset", 
        dtype_id = "H5T_STD_I64LE", 
        h5space = rhdf5::H5Screate_simple(length(indices$bin1_offset))
    )
    rhdf5::H5Dwrite(
        did, 
        as.integer(indices$bin1_offset)
    )
    rhdf5::H5Dclose(did)
    #  chrom_offset 
    did <- rhdf5::H5Dcreate(
        h5loc = fid, 
        name = "indexes/chrom_offset", 
        dtype_id = "H5T_STD_I64LE", 
        h5space = H5Screate_simple(length(indices$chrom_offset))
    )
    rhdf5::H5Dwrite(
        did, 
        as.integer(indices$chrom_offset)
    )
    rhdf5::H5Dclose(did)
    rhdf5::H5Gclose(gid)

    ## Close all handles
    rhdf5::h5closeAll()

    return(output)
}
