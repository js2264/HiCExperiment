% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{HiCExperiment utils}
\alias{HiCExperiment utils}
\alias{splitCoords}
\alias{coords2char}
\alias{char2coords}
\alias{sortPairs}
\alias{asGInteractions}
\alias{sdiag}
\alias{.df_to_symmmat}
\alias{distance_decay}
\title{Utils functions}
\usage{
splitCoords(coords)

coords2char(coords, big.mark = ",")

char2coords(char)

sortPairs(pairs)

asGInteractions(df)

sdiag(A, k = 0)

.df_to_symmmat(diag, score)

distance_decay(dump, threshold = 1e+10)
}
\arguments{
\item{coords}{coords}

\item{big.mark}{big.mark}

\item{char}{char (e.g. "II:30000-50000" or "II:30000-50000|II:60000-80000")}

\item{pairs}{pairs}

\item{df}{df}

\item{A}{matrix}

\item{k}{secondary diagonal k}

\item{diag}{vector of distances to diagonal}

\item{score}{scores to parse into symmetrical matrix}

\item{dump}{dumped contacts, e.g. from .dumpCool}

\item{threshold}{maximum distance to compute distance decay for}
}
\value{
Reformatted coordinates or GInteractions
}
\description{
Utilities to facilitate parsing/handling of coordinates, GInteractions,
Pairs, ...
}
