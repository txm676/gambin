#' @name moths
#' @aliases moths
#' @title Williams' Rothamsted moth data
#' @description Macro-Lepidoptera captured in a light trap at Rothamsted Experimental Station during 1935.
#' @docType data
#' @format A numerical vector with the abundance of 195 moth species.
#' @source Williams, C.B. (1964) Patterns in the balance of nature. Academic Press, London.
#' @examples
#' data(moths, package = "gambin")
NULL


#' @name fly
#' @aliases fly
#' @title Brazilian Horse Fly Data
#' @description Horse flies captured using various sampling methods at different sites
#'   across Brazil.
#' @docType data
#' @format A list with two elements. The first element contains a numerical
#'   vector with the abundance of 164 fly species sampled at various sites
#'   across Brazil. The second element contains a numerical vector with the
#'   abundance of 58 fly specie sampled at a single site within Brazil using
#'   just canopy traps.
#' @source This package.
#' @examples
#' data(fly, package = "gambin")
NULL


#' @name categ
#' @aliases categ
#' @title Simulated bird SAD dataset with species classification data
#' @description A randomly generated bird SAD dataset where each species has
#'   been randomly classified according to its origin (native, exotic or
#'   invasive).
#' @docType data
#' @format A dataframe with three columns: 1) 'abundances' = the abundance of
#'   each species, 2) 'species' = the species names, and 3) 'status' the species
#'   origin classification. In regards to (3) each species is classified as
#'   either native (N), exotic (E) or invasive (I).
#' @source This package.
#' @examples
#' data(categ, package = "gambin")
NULL