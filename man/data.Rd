% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data}
\alias{data}
\alias{centros_yeast}
\alias{contacts_yeast}
\alias{contacts_yeast_eco1}
\title{Example datasets provided in \code{HiCExperiment} & \code{HiContactsData}}
\format{
An object of class \code{"GRanges"}.
}
\source{
HiContacts
}
\usage{
data(centros_yeast)

contacts_yeast(full = FALSE)

contacts_yeast_eco1(full = FALSE)
}
\arguments{
\item{full}{Whether to import all interactions}
}
\description{
Example datasets provided in \code{HiCExperiment} & \code{HiContactsData}
}
\examples{
data(centros_yeast)
centros_yeast
contacts_yeast()
}
