% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import-methods.R
\name{import-methods}
\alias{import-methods}
\alias{import}
\alias{import,CoolFile-method}
\alias{import,HicFile-method}
\alias{import,HicproFile-method}
\alias{import,PairsFile-method}
\alias{import,CoolFile,ANY,ANY-method}
\alias{import,HicFile,ANY,ANY-method}
\alias{import,HicproFile,ANY,ANY-method}
\alias{import,PairsFile,ANY,ANY-method}
\title{HiCExperiment import methods}
\usage{
import(con, format, text, ...)
}
\arguments{
\item{con}{Path or connection to a cool, mcool, .hic or HiC-Pro derived files.
Can also be a path to a pairs file.}

\item{format}{The format of the output. If missing and 'con' is a filename,
the format is derived from the file extension.
This argument is unnecessary when files are directly provided as
\code{CoolFile}, \code{HicFile}, \code{HicproFile} or \code{PairsFile}.}

\item{text}{If 'con' is missing, this can be a character vector directly
providing the string data to import.}

\item{...}{e.g. \code{resolution = ...}; parameters to pass to
format-specific methods.}
}
\value{
A \code{HiCExperiment} or \code{GInteractions} object
}
\description{
Import methods for data structures implemented in the HiCExperiment package.
HiCExperiment package implements methods to faciliate
the import of Hi-C files (.(m)cool, .hic, HiC-Pro derived)
and pairs files in R, as HiCExperiment or GInteractions objects.
}
\examples{
# ---- Importing .(m)cool files 
mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
import(mcool_path, resolution = 16000, focus = 'XVI', format = 'cool')

# ---- Importing .hic files 
hic_path <- HiContactsData::HiContactsData('yeast_wt', 'hic')
import(hic_path, resolution = 16000, focus = 'XVI', format = 'hic')

# ---- Importing HiC-Pro files 
#hicpro_matrix_path <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_matrix')
#hicpro_bed_path <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_bed')
#import(hicpro_matrix_path, bed = hicpro_bed_path, format = 'hicpro')
}
