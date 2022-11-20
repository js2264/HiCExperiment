#' @title Generic functions
#' 
#' @name AllGenerics
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
#' 
#' @description
#' 
#' Generics functions created in HiCExperiment package. 
#' 
#' @param x Passed to corresponding method
#' @param name Passed to corresponding method
#' @param value Passed to corresponding method
NULL

#' @export 

setGeneric("resolutions", function(x) {standardGeneric("resolutions")})

#' @export 

setGeneric("resolution", function(x) {standardGeneric("resolution")})

#' @export 

setGeneric("focus", function(x) {standardGeneric("focus")})

#' @export 

setGeneric("focus<-", function(x, value) {standardGeneric("focus<-")})

#' @export 

setGeneric("scores", function(x, name) {standardGeneric("scores")})

#' @export 

setGeneric("scores<-", function(x, name, value) {standardGeneric("scores<-")})

#' @export 

setGeneric("topologicalFeatures", function(x, name) {standardGeneric("topologicalFeatures")})

#' @export 

setGeneric("topologicalFeatures<-", function(x, name, value) {standardGeneric("topologicalFeatures<-")})

#' @export 

setGeneric("pairsFile", function(x, name) {standardGeneric("pairsFile")})

#' @export 

setGeneric("pairsFile<-", function(x, value) {standardGeneric("pairsFile<-")})

#' @export 

setGeneric("metadata<-", function(x, value) {standardGeneric("metadata<-")})

#' @export 

setGeneric("bins", function(x) {standardGeneric("bins")})
