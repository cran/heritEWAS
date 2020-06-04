#' Calculate carrier probabilities for the most heritable methylation sites
#'
#' For each person in `dat` and each methylation site in `top_probes`,
#' this function calculates the probability that the person carries a rare
#' mutation at a hypothetical genetic locus that affects methylation
#' at the methylation site.
#'
#' @inheritParams genotype_combinations
#' @inheritParams ML_estimates
#' @param top_probes  A data frame, usually the output of \code{\link{ML_estimates}} restricted to
#' the most heritable methylation sites (those with the highest values of \eqn{\Delta}l).
#' See \code{\link{ML_estimates}} and the example below for more details.
#'
#' @param ncores The number of cores to be used, with `ncores = 1` (the
#' default) corresponding to non-parallel computing.  When `ncores > 1`,
#' the `parallel` package is used to parallelize the calculation.
#'
#' @param M_values A matrix of M-values, with rows corresponding to
#' methylation sites and columns corresponding to people.
#'
#' @return A data frame containing the carrier probabilities described in (Joo et al., 2018),
#' with rows of the data frame corresponding to the people in `dat` and columns corresponding to
#' the methylation sites in `top_probes`.  This calculation is based on the Mendelian model
#' of (Joo et al., 2018) with parameter values taken from `top_probes`.
#'
#' @export
#'
#' @examples
#' str(ped)
#' str(M_values)
#'
#' # Calculate genotype probabilities
#' typed_genos <- genotype_combinations(ped)
#' str(typed_genos)
#'
#' \donttest{
#' # Compute Delta l
#' MLEs <- ML_estimates(typed_genos, M_values, ncores = 4)
#'
#' # Select top probes
#' top_probes <- MLEs[MLEs$delta.l > 10, ]
#'
#' # Calculate carrier probabilities
#' CP <- carrier_probabilities(ped, M_values, top_probes, ncores = 4)
#' str(CP)
#' }
#'
#' @references
#' Joo JE, Dowty JG, Milne RL, Wong EM, Dugué PA, English D, Hopper JL,
#' Goldgar DE, Giles GG, Southey MC, kConFab.  Heritable DNA methylation marks
#' associated with susceptibility to breast cancer.  Nat Commun. 2018
#' Feb 28;9(1):867. \url{https://doi.org/10.1038/s41467-018-03058-6}
#'
carrier_probabilities <- function(dat, M_values, top_probes, ncores = 1) {

  # Check data
  check_dat(dat)

  # Read M-values
  M <- M_values

  # Restrict M to top probes
  w <- which(rownames(M) %in% top_probes$methylation.site)
  M <- M[w, ]

  # Make a copy of the original ID
  dat$ID.orig <- dat$indiv

  # Covert ID
  dat <- convert_IDs(dat, convert.IDs.numeric = TRUE)
  datp <- dat

  # Find core family
  ufID <- unique(dat$family)
  core.fams <- vector("list", length(ufID))
  for (i in seq_along(ufID)) {
    fam <- dat[dat$family == ufID[i],
               names(dat) %in% c("indiv", "mother", "father",
                                 "ID.orig", "typed")]
    core <- trim_pedigree(fam)
    core.fams[[i]] <- core
  }
  names(core.fams) <- ufID

  # Find the columns of M corresponding to the columns of each element of
  # core.fams
  find_indices_core <- function(fam) {
    x <- fam$ID.orig
    index <- numeric(length(x))
    for (i in seq_along(x)) {
      if (!any(colnames(M) == x[i])) {
        index[i] <- NA
      } else {
        index[i] <- which(colnames(M) == x[i])[1]
      }
    }
    index
  }
  indices.core.fams <- lapply(core.fams, find_indices_core)

  # Find the rows of the full family corresponding to each member of the core
  # families
  find_indices_map_from_core <- function(fam.name) {
    fam <- core.fams[names(core.fams) == fam.name][[1]]
    x <- fam$ID.orig
    y <- datp$ID.orig[datp$family == fam.name]
    index <- numeric(length(x))
    for (i in 1:length(x)) {
      if (!any(y == x[i])) {
        index[i] <- NA
      } else {
        index[i] <- which(y == x[i])[1]
      }
    }
    index
  }
  indices.map.from.core.fams <- lapply(names(core.fams),
                                       find_indices_map_from_core)

  # Find the rows of datp corresponding to each member of the full families
  find_indices_map_to_datp <- function(fam.name) {
    full <- datp[datp$family == fam.name, ]
    x <- full$ID.orig
    y <- datp$ID.orig
    rightfam <- (datp$family == fam.name)
    index <- numeric(length(x))
    for (i in 1:length(x)) {
      if (!any(y == x[i] & rightfam)) {
        index[i] <- NA
      } else {
        index[i] <- which(y == x[i] & rightfam)[1]
      }
    }
    index
  }
  indices.map.to.datp <- lapply(names(core.fams), find_indices_map_to_datp)

  # Compute genotype combinations for the core families
  core.genos <- calc_core_probs(core.fams)

  # Read top probes
  top <- top_probes

  # Calculate the carrier probabilities corresponding to a given probe
  # Loop over the probes

  calc_for_one_probe <- function(j.for.probe) {

    # Set the parameters to the ML estimates
    mu0 <- top$mu0.mendel[j.for.probe]
    sd0 <- top$sd0.mendel[j.for.probe]
    mu1 <- top$mu1.mendel[j.for.probe]
    sd1 <- top$sd1.mendel[j.for.probe]
    c(mu0, sd0, mu1, sd1)

    # Loop over the families
    prob.probej <- rep(NA, nrow(datp))
    for (i.for.fam in seq_along(core.fams)) {

      # Initialise
      fam <- core.fams[[i.for.fam]]
      curr.fam.name <- names(core.fams)[i.for.fam]

      # Read in M-values
      ind <- indices.core.fams[[i.for.fam]]
      mval <- M[rownames(M) == top$methylation.site[j.for.probe], ind]

      # All possible genotype combinations for the core family
      genos <- core.genos[[i.for.fam]]

      # Calculate the densities corresponding to the M-values
      f0 <- function(x) dnorm(x, mu0, sd0)
      f1 <- function(x) dnorm(x, mu1, sd1)
      y <- c(f0(mval), f1(mval))
      p01 <- matrix(y, length(mval), 2)
      p01[is.na(p01)] <- 1

      # Calculate the probability for each genotype combination
      p <- genos$p
      for (k in 1:nrow(fam)) {
        indices <- genos[, k] + 1
        p <- p * p01[k, indices]
      }

      # Calculate the carrier probabilities for each person in the core family
      carr.prob <- rep(-1, nrow(fam))
      for (k in 1:nrow(fam)) {
        p0 <- sum(p[genos[, k] == 0])
        p1 <- sum(p[genos[, k] == 1])
        carr.prob[k] <- p1 / (p0 + p1)
      }
      p01[fam$typed == 1, ]
      carr.prob[fam$typed == 1]

      # Extend the carrier probabilities to people outside of the core family
      full <- datp[datp$family == curr.fam.name, ]
      cp <- rep(NA, nrow(full))
      ind <- indices.map.from.core.fams[[i.for.fam]]
      cp[ind] <- carr.prob
      # Assume founders not ancestoprs of some core person have no chance of
      # carrying a mutation
      corefamIDs <- full$indiv[full$ID.orig %in% fam$ID.orig]
      poss.carriers <- descendents(ancestors2(corefamIDs, full), full)
      cp[!(full$indiv %in% poss.carriers)] <- 0
      data.frame(full, cp, core = as.numeric(full$indiv %in% corefamIDs))

      # A function to extend the people with non-missing carrier probabilities
      expand <- function(cp, full) {
        for (i in 1:nrow(full)) {
          if (is.na(cp[i])) {
            children <- full$mother == full$indiv[i] |
              full$father == full$indiv[i]
            children <- children & !is.na(cp)
            if (any(children)) {
              if (sum(children) >= 2) {
                warning("Filling in carrier probs with two or more children of",
                        " ", full$indiv[i], call. = TRUE)
              }
              cp[i] <- max(cp[children]) / 2
            }
            mum <- full$indiv == full$mother[i]
            dad <- full$indiv == full$father[i]
            if (any(mum) & any(dad)) {
              if (!is.na(cp[mum]) & !is.na(cp[dad])) {
                cp[i] <- (cp[mum] + cp[dad]) / 2
              }
            }
          }
        }
        cp
      }

      # Calculate carrier probabilities for everyone in the family
      n.old <- sum(is.na(cp)) + 1
      while (n.old > sum(is.na(cp))) {
        n.old <- sum(is.na(cp))
        cp <- expand(cp, full)
      }

      # Add them to the probabilities for all families
      ind <- indices.map.to.datp[[i.for.fam]]
      prob.probej[ind] <- round(cp, 4)
    }
    prob.probej
  }

  cl <- parallel::makeCluster(ncores)
  export_list <- list("top", "core.fams", "datp", "M",
                      "indices.map.from.core.fams",
                      "indices.core.fams", "core.genos",
                      "descendents", "ancestors2",
                      "indices.map.to.datp")

  parallel::clusterExport(cl, export_list, envir = environment())
  on.exit(parallel::stopCluster(cl), add = TRUE)
  probe_probs <- parallel::parSapply(cl, 1:nrow(top), calc_for_one_probe)
  colnames(probe_probs) <- top$methylation.site

  data.frame(dat, probe_probs)

}
