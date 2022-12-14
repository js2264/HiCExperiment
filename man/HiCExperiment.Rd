% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HiCExperiment-class.R,
%   R/HiCExperiment-methods.R
\docType{class}
\name{HiCExperiment}
\alias{HiCExperiment}
\alias{HiCExperiment-class}
\alias{.HiCExperimentFromCoolFile}
\alias{.HiCExperimentFromHicFile}
\alias{.HiCExperimentFromHicproFile}
\alias{resolutions,HiCExperiment-method}
\alias{resolution,HiCExperiment-method}
\alias{focus,HiCExperiment-method}
\alias{focus<-,HiCExperiment-method}
\alias{focus<-,HiCExperiment,character-method}
\alias{zoom,HiCExperiment,numeric-method}
\alias{refocus,HiCExperiment,character-method}
\alias{scores,HiCExperiment-method}
\alias{scores,HiCExperiment,missing-method}
\alias{scores,HiCExperiment,character-method}
\alias{scores,HiCExperiment,numeric-method}
\alias{scores<-,HiCExperiment-method}
\alias{scores<-,HiCExperiment,character,numeric-method}
\alias{topologicalFeatures,HiCExperiment-method}
\alias{topologicalFeatures,HiCExperiment,missing-method}
\alias{topologicalFeatures,HiCExperiment,character-method}
\alias{topologicalFeatures,HiCExperiment,numeric-method}
\alias{topologicalFeatures<-,HiCExperiment-method}
\alias{topologicalFeatures<-,HiCExperiment,character,GRangesOrGInteractions-method}
\alias{pairsFile,HiCExperiment-method}
\alias{pairsFile<-,HiCExperiment-method}
\alias{pairsFile<-,HiCExperiment,character-method}
\alias{metadata<-,HiCExperiment-method}
\alias{metadata<-,HiCExperiment,list-method}
\alias{fileName,HiCExperiment-method}
\alias{interactions,HiCExperiment-method}
\alias{interactions<-,HiCExperiment-method}
\alias{length,HiCExperiment-method}
\alias{[,HiCExperiment-method}
\alias{[,HiCExperiment,numeric,ANY,ANY-method}
\alias{[,HiCExperiment,logical,ANY,ANY-method}
\alias{[,HiCExperiment,character,ANY,ANY-method}
\alias{seqinfo,HiCExperiment-method}
\alias{bins,HiCExperiment-method}
\alias{anchors,HiCExperiment-method}
\alias{regions,HiCExperiment-method}
\alias{show,HiCExperiment-method}
\alias{coerce,HiCExperiment,GInteractions-method}
\alias{coerce,HiCExperiment,ContactMatrix-method}
\alias{coerce,HiCExperiment,matrix-method}
\alias{coerce,HiCExperiment,data.frame-method}
\alias{as.matrix,HiCExperiment-method}
\alias{as.data.frame,HiCExperiment-method}
\alias{interactions<-,HiCExperiment,GInteractions-method}
\alias{as}
\title{\code{HiCExperiment} S4 class}
\usage{
HiCExperiment(
  file,
  resolution = NULL,
  focus = NULL,
  metadata = list(),
  topologicalFeatures = S4Vectors::SimpleList(compartments = GenomicRanges::GRanges(),
    borders = GenomicRanges::GRanges(), loops =
    InteractionSet::GInteractions(GenomicRanges::GRanges(), GenomicRanges::GRanges()),
    viewpoints = GenomicRanges::GRanges()),
  pairsFile = NULL,
  ...
)

.HiCExperimentFromCoolFile(
  file,
  resolution = NULL,
  focus = NULL,
  metadata = list(),
  topologicalFeatures = S4Vectors::SimpleList(compartments = GenomicRanges::GRanges(),
    borders = GenomicRanges::GRanges(), loops =
    InteractionSet::GInteractions(GenomicRanges::GRanges(), GenomicRanges::GRanges()),
    viewpoints = GenomicRanges::GRanges()),
  pairsFile = NULL
)

.HiCExperimentFromHicFile(
  file,
  resolution = NULL,
  focus = NULL,
  metadata = list(),
  topologicalFeatures = S4Vectors::SimpleList(compartments = GenomicRanges::GRanges(),
    borders = GenomicRanges::GRanges(), loops =
    InteractionSet::GInteractions(GenomicRanges::GRanges(), GenomicRanges::GRanges()),
    viewpoints = GenomicRanges::GRanges()),
  pairsFile = NULL
)

.HiCExperimentFromHicproFile(
  file,
  bed,
  metadata = list(),
  topologicalFeatures = S4Vectors::SimpleList(compartments = GenomicRanges::GRanges(),
    borders = GenomicRanges::GRanges(), loops =
    InteractionSet::GInteractions(GenomicRanges::GRanges(), GenomicRanges::GRanges()),
    viewpoints = GenomicRanges::GRanges()),
  pairsFile = NULL
)

\S4method{resolutions}{HiCExperiment}(x)

\S4method{resolution}{HiCExperiment}(x)

\S4method{focus}{HiCExperiment}(x)

\S4method{focus}{HiCExperiment,character}(x) <- value

\S4method{zoom}{HiCExperiment,numeric}(x, resolution)

\S4method{refocus}{HiCExperiment,character}(x, focus)

\S4method{scores}{HiCExperiment,missing}(x)

\S4method{scores}{HiCExperiment,character}(x, name)

\S4method{scores}{HiCExperiment,numeric}(x, name)

\S4method{scores}{HiCExperiment,character,numeric}(x, name) <- value

\S4method{topologicalFeatures}{HiCExperiment,missing}(x)

\S4method{topologicalFeatures}{HiCExperiment,character}(x, name)

\S4method{topologicalFeatures}{HiCExperiment,numeric}(x, name)

\S4method{topologicalFeatures}{HiCExperiment,character,GRangesOrGInteractions}(x, name) <- value

\S4method{pairsFile}{HiCExperiment}(x)

\S4method{pairsFile}{HiCExperiment,character}(x) <- value

\S4method{metadata}{HiCExperiment,list}(x) <- value

\S4method{fileName}{HiCExperiment}(object)

\S4method{interactions}{HiCExperiment}(x)

\S4method{interactions}{HiCExperiment,GInteractions}(x) <- value

\S4method{length}{HiCExperiment}(x)

\S4method{[}{HiCExperiment,numeric,ANY,ANY}(x, i)

\S4method{[}{HiCExperiment,logical,ANY,ANY}(x, i)

\S4method{[}{HiCExperiment,character,ANY,ANY}(x, i)

\S4method{seqinfo}{HiCExperiment}(x)

\S4method{bins}{HiCExperiment}(x)

\S4method{anchors}{HiCExperiment}(x)

\S4method{regions}{HiCExperiment}(x)

\S4method{show}{HiCExperiment}(object)

\S4method{as.matrix}{HiCExperiment}(x)

\S4method{as.data.frame}{HiCExperiment}(x)
}
\arguments{
\item{file}{CoolFile or plain path to a Hi-C contact file}

\item{resolution}{Resolution to use with the Hi-C contact file}

\item{focus}{Chromosome coordinates for which
interaction counts are extracted from the Hi-C contact file, provided
as a character string (e.g. "II:4000-5000"). If not provided,
the entire Hi-C contact file will be imported.}

\item{metadata}{list of metadata}

\item{topologicalFeatures}{topologicalFeatures provided as a named SimpleList}

\item{pairsFile}{Path to an associated .pairs file}

\item{...}{Extra arguments}

\item{bed}{Path to regions file generated by HiC-Pro}

\item{x}{A \code{HiCExperiment} object.}

\item{value}{...}

\item{name}{...}

\item{object}{A \code{HiCExperiment} object.}

\item{i}{a range or boolean vector.}
}
\value{
An \code{HiCExperiment} object.
}
\description{
The \code{HiCExperiment} class describes Hi-C contact files imported in R, either
through the \code{HiCExperiment} constructor function or using the \code{import}
method implemented by \code{HiCExperiment} package.

HiCExperiment methods.
}
\section{Slots}{

\describe{
\item{\code{fileName}}{Path of Hi-C contact file}

\item{\code{focus}}{Chr. coordinates for which interaction counts are extracted
from the Hi-C contact file.}

\item{\code{resolutions}}{Resolutions available in the Hi-C contact file.}

\item{\code{resolution}}{Current resolution}

\item{\code{interactions}}{Genomic Interactions extracted from the Hi-C contact file}

\item{\code{scores}}{Available interaction scores.}

\item{\code{topologicalFeatures}}{Topological features associated with the dataset
(e.g. loops (\<Pairs\>), borders (\<GRanges\>),
viewpoints (\<GRanges\>), etc...)}

\item{\code{pairsFile}}{Path to the .pairs file associated with the Hi-C contact file}

\item{\code{metadata}}{metadata associated with the Hi-C contact file.}
}}

\examples{
mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
contacts_yeast <- HiCExperiment(mcool_path, resolution = 16000, focus = 'II')
contacts_yeast
resolutions(contacts_yeast)
resolution(contacts_yeast)
focus(contacts_yeast)
scores(contacts_yeast)
tail(scores(contacts_yeast, 1))
tail(scores(contacts_yeast, 'balanced'))
scores(contacts_yeast, 'test') <- runif(length(contacts_yeast))
tail(scores(contacts_yeast, 'test'))
data(centros_yeast)
topologicalFeatures(contacts_yeast, 'centromeres') <- centros_yeast
topologicalFeatures(contacts_yeast, 1)
topologicalFeatures(contacts_yeast, 'centromeres')
pairsFile(contacts_yeast)
fileName(contacts_yeast)
interactions(contacts_yeast)
length(contacts_yeast)
contacts_yeast[seq_len(10)]
seqinfo(contacts_yeast)
bins(contacts_yeast)
anchors(contacts_yeast)
regions(contacts_yeast)
contacts_yeast
as(contacts_yeast, 'GInteractions')
as(contacts_yeast, 'ContactMatrix')
as(contacts_yeast, 'matrix')[seq_len(10), seq_len(10)]
as(contacts_yeast, 'data.frame')[seq_len(10), seq_len(10)]
}
