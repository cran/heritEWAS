#' Calculate genotype probabilities
#'
#' This function computes the joint probabilities of the possible genotypes of
#' selected family members within each family, as an intermediate calculation
#' for use in \code{\link{ML_estimates}}.
#'
#' @param dat A data frame with rows corresponding to people and columns
#' corresponding to (at least) the following variables, which will be coerced to
#' `character` type:
#'
#' * `family` (family ID), an identifier for each person's family, constant
#'    within families
#' * `indiv` (individual ID), an identifier for each person, with no duplicates
#'    across the dataset
#' * `mother` (mother ID), the individual ID of each person's mother, or missing
#'    (`NA`) for founders
#' * `father` (father ID), the individual ID of each person's father, or missing
#'    (`NA`) for founders
#' * `typed` (epi-genotyped), equal to `1` for people with methylation data and
#'    `0` for all others
#'
#' @param ncores The number of cores to be used, with `ncores = 1` (the
#' default) corresponding to non-parallel computing.  When `ncores > 1`,
#' the `parallel` package is used to parallelize the calculation of the joint
#' probabilities of the possible genotype combinations.
#'
#' @param verbose \code{FALSE} if user wants to suppress messages. Default is
#' \code{TRUE}.
#'
#' @details
#' Currently, there is a maximum of 20 people per family with `typed = 1`,
#' due to the need to store all genotype combinations for the typed people.
#' It is possible to break families with more than 20 typed people into
#' smaller families, though this is not ideal.
#'
#' Each family within `dat` should be a complete pedigree, meaning that each
#' (non-missing) mother or father ID should correspond to a row, and each person
#' should either have both parent IDs missing (if a founder) or non-missing
#' (if a non-founder).  No family should contain a pedigree loop, such as those
#' caused by inbreeding or by two sisters having children with two brothers
#' from an unrelated family.
#'
#' @return A named list, with names equal to the different family IDs and with
#' each element of the list being a data frame specifying the possible genotypes
#' of selected family members (those with `dat$typed = 1`) within each family,
#' and the joint probability of each genotype combination.  For each individual,
#' the possible genotypes are `0`, corresponding to the wildtype, and `1`,
#' for carriers of a rare genetic variant at an autosomal locus. The calculation
#' makes the same assumptions as in (Joo et al., 2018), including that the
#' genetic variant is so rare that at most one founder is a carrier.
#'
#' @export
#'
#' @examples
#' # Load family data
#' data(ped)
#'
#' # Calculate genotype probabilites
#' typed_genos <- genotype_combinations(ped)
#' str(typed_genos)
#'
#' @references
#' Joo JE, Dowty JG, Milne RL, Wong EM, Dugué PA, English D, Hopper JL,
#' Goldgar DE, Giles GG, Southey MC, kConFab.  Heritable DNA methylation marks
#' associated with susceptibility to breast cancer.  Nat Commun. 2018
#' Feb 28;9(1):867. \url{https://doi.org/10.1038/s41467-018-03058-6}
#'
genotype_combinations <- function(dat, ncores = 1, verbose = TRUE) {

  check_dat(dat)
  dat$family <- as.character(dat$family)
  dat$indiv <- as.character(dat$indiv)
  dat$mother <- as.character(dat$mother)
  dat$father <- as.character(dat$father)

  dat$mother[which(is.na(dat$mother))] <- ""
  dat$father[which(is.na(dat$father))] <- ""

  ids <- unique(dat$family)

  if (any(!ids %in% unique(dat$family))) {
    stop("family cannot be found in the pedigree data.")
  }

  tg_list <- vector("list", length(ids))

  j <- 1
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  for (i in ids) {

    if (verbose) {
      cat("Family ", j, " (", i, ")", ": ", sep = "")
    }

    fam <- dat[dat$family == i,
               names(dat) %in% c("indiv", "mother", "father", "ID.orig",
                                 "typed")]
    if (sum(fam$typed) > 20) {
      stop("Maximum of 20 people per family with typed = 1 is allowed, consider splitting pedigrees")
    }

    tg_list[[j]] <- genotype_combinations_f(fam, ncores, cl, verbose)
    j <- j + 1
  }

  names(tg_list) <- ids

  # Adjusting genotype probabilites for Mendelian model
  tg_list <- adjust_geno_probs(tg_list)

  tg_list

}
