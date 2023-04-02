#' @title Generic functions
#' 
#' @name AllGenerics
#' @aliases availableResolutions
#' @aliases availableChromosomes
#' @aliases resolutions
#' @aliases resolution
#' @aliases focus
#' @aliases focus<-
#' @aliases scores
#' @aliases scores<-
#' @aliases topologicalFeatures
#' @aliases topologicalFeatures<-
#' @aliases pairsFile
#' @aliases pairsFile<-
#' @aliases metadata<-
#' @aliases bins
#' @aliases slices
#' @aliases zoom
#' @aliases refocus
#' @aliases cis
#' @aliases trans
#' @aliases bin
#' 
#' @description
#' 
#' Generics functions created in HiCExperiment package. 
#' 
#' @param x Passed to corresponding method
#' @param name Passed to corresponding method
#' @param value Passed to corresponding method
#' @param ... Passed to corresponding method
#' @keywords internal
NULL

setGeneric("availableResolutions", function(x, ...) {standardGeneric("availableResolutions")})
setGeneric("availableChromosomes", function(x, ...) {standardGeneric("availableChromosomes")})
setGeneric("resolutions", function(x) {standardGeneric("resolutions")})
setGeneric("resolution", function(x) {standardGeneric("resolution")})
setGeneric("focus", function(x) {standardGeneric("focus")})
setGeneric("focus<-", function(x, value) {standardGeneric("focus<-")})
setGeneric("scores", function(x, name) {standardGeneric("scores")})
setGeneric("scores<-", function(x, name, value) {standardGeneric("scores<-")})
setGeneric("topologicalFeatures", function(x, name) {standardGeneric("topologicalFeatures")})
setGeneric("topologicalFeatures<-", function(x, name, value) {standardGeneric("topologicalFeatures<-")})
setGeneric("pairsFile", function(x, name) {standardGeneric("pairsFile")})
setGeneric("pairsFile<-", function(x, value) {standardGeneric("pairsFile<-")})
setGeneric("metadata<-", function(x, value) {standardGeneric("metadata<-")})
setGeneric("bins", function(x) {standardGeneric("bins")})
setGeneric("slices", function(x, name) {standardGeneric("slices")})
setGeneric("zoom", function(x, resolution) {standardGeneric("zoom")})
setGeneric("refocus", function(x, focus) {standardGeneric("refocus")})
setGeneric("cis", function(x) {standardGeneric("cis")})
setGeneric("trans", function(x) {standardGeneric("trans")})
setGeneric("bin", function(x, ...) {standardGeneric("bin")})
