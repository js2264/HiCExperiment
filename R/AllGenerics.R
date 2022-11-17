#' @title Generic functions
#' 
#' @name AllGenerics
#' 
#' @description
#' 
#' Generics functions created in HiCExperiment package. 
#' 
#' @param x Passed to corresponding method
#' @param name Passed to corresponding method
#' @param value Passed to corresponding method
NULL

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
