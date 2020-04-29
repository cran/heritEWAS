#' heritEWAS: Identify heritable DNA methylation marks from families
#'
#' `heritEWAS` implements the statistical method of (Joo et al., 2018)
#' to efficiently scan the genome for DNA methylation marks that are heritable.
#' This method looks for Mendelian patterns of inheritance at methylation sites,
#' and is based only on the relationship structure of (large) families and on
#' the methylation data of some individuals within these families.
#' In particular, genotype data is not required.
#'
#' The methylation data will typically come from
#' an epigenome-wide association study (EWAS), and `heritEWAS` is optimised for
#' such data, e.g. the most time-consuming part of the calculation is performed
#' once, then the output of this calculation is re-used for all methylation
#' sites. However, the code can also be run on methylation data at a single
#' genomic location.
#'
#' For each methylation site, the code computes a measure of heritability called
#' Delta l. It is not clear what relationship exists between this measure
#' and classical measures of heritability, such as those estimated by twin
#' studies.
#'
#' The methods implemented in this package are described briefly in
#' \code{\link{compute_deltal}} and fully in (Joo et al., 2018).
#'

#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

