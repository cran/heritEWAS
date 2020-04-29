calc_geno_probs_f <- function(fam, ncores = 1, cl = NULL, verbose) {

  # Make a copy of the original ID
  fam$ID.orig <- fam$indiv

  # Trim the pedigrees w.r.t typed people
  fam <- trim_pedigree(fam, reorder = TRUE)

  # Check if only one typed person in the family
  if (nrow(fam) == 1) {
    out <- data.frame(p = c(1, 1), x = c(0, 1))
    names(out)[2] <- fam$ID.orig[1]
    return(out)
  }

  # Find the founders, typed people and their intersection
  founders <- fam$indiv[fam$mother == 0 & fam$father == 0]
  typed <- fam$indiv[fam$typed == 1]
  common <- intersect(founders, typed)

  # Possible genotypes for typed people
  typed.geno <- generate_typed_geno(fam)

  if (verbose) {
    cat("Number of genotype combinations to check:", nrow(typed.geno), "\n")
  }

  # Possible genotypes for founders to order 1
  founder.geno.deg1 <- diag(length(founders))

  # Possible genotype for remaining people who are neither founders nor typed
  remaining <- setdiff(1:nrow(fam), c(typed, founders))
  if (length(remaining) == 0) {
    remaining.geno <- matrix(0, 0, 0)
  } else {
    remaining.geno <- geno.comb(0:1, length(remaining))
  }

  # Create matrix of transmission probablities
  trans.matrix <- create_trans_matrix()

  # Calculate the transmission probabilities
  trans.prob <- function(g, trans.matrix, fam) {
    # Given a vector g = (g1, ..., gn) of unphased genotypes,
    # one for each member of fam with each gi = 0, 1 or 2,
    # calculate the transmission probabilities, i.e. P(g)/Prior(g).
    # This function assumes the people in fam are numbered consecutively
    # and either the mum and dad IDs are both 0 or neither are.

    # Reorder the people in fam consecutively
    # fam <- convert.IDs(cbind(family = 1, fam),
    #                    convert.IDs.numeric = TRUE)[, -1]
    gm <- g[fam$mother]
    gd <- g[fam$father]
    go <- g[fam$mother != 0 & fam$father != 0]
    y <- 1 + 9 * go + 3 * gm + gd
    p <- trans.matrix$p[y]
    prod(p)
  }

  # Calculate the probability for each genotype combinaton
  p.trans.deg1 <- function(gtyped) {
    gtyped <- as.numeric(gtyped)
    g <- rep(-1, nrow(fam))
    g[typed] <- gtyped
    trans <- 0

    # If the founders and typed people overlap then restrict founder.geno.deg1
    # so it is consistent with gtyped
    fgd1 <- founder.geno.deg1

    if (length(common) > 0) {
      keep <- rep(TRUE, nrow(fgd1))
      for (i in 1:length(common)) {
        j <- which(founders == common[i])[1]
        k <- which(typed == common[i])[1]
        keep <- keep & (fgd1[, j] == gtyped[k])
      }
      if (sum(keep) == 0) {
        return(0)
      }
      fgd1 <- fgd1[keep, ]

      if (sum(keep) == 1) fgd1 <- t(as.matrix(as.vector(fgd1)))
    }

    if (length(remaining) == 0) {
      for (i in 1:nrow(fgd1)) {
        g[founders] <- fgd1[i, ]
        trans <- trans + trans.prob(g, trans.matrix, fam)
      }
    } else if (nrow(remaining.geno) > 0) {
      for (i in 1:nrow(fgd1)) {
        for (j in 1:nrow(remaining.geno)) {
          g[founders] <- fgd1[i, ]
          g[remaining] <- as.numeric(remaining.geno[j, ])
          trans <- trans + trans.prob(g, trans.matrix, fam)
        }
      }
    }
    # Prior(g) = 2 * p^k * (1-p)^(2*f-k) when k==1
    trans <- 2 * trans
    trans
  }

  if (ncores == 1) {
    prob <- apply(typed.geno, 1, p.trans.deg1)
  } else {
    #cl <- parallel::makeCluster(ncores)
    #on.exit(parallel::stopCluster(cl), add = TRUE)
    prob <- parallel::parApply(cl, typed.geno, 1, p.trans.deg1)
  }

  # Output file
  out <- data.frame(prob[prob != 0], typed.geno[prob != 0, ])
  names(out) <- c("p", fam$ID.orig[typed])
  row.names(out) <- NULL

  out

}

