calc_core_probs <- function(core.fams) {

  # # Make a copy of original ID and convert ID into 1 to n
  # dat$ID.orig <- dat$indiv
  # dat <- convert_IDs(dat, convert.IDs.numeric = TRUE)
  #
  # # Find core family
  # ufID <- unique(dat$family)
  # core.fams <- vector("list", length(ufID))
  # for (i in seq_along(ufID)) {
  #   fam <- dat[dat$family == ufID[i],
  #               names(dat) %in% c("indiv", "mother", "father",
  #                                 "ID.orig", "typed")]
  #   core <- trim_pedigree(fam)
  #   core.fams[[i]] <- core
  # }
  # names(core.fams) <- ufID

  # Loop over the families
  results <- vector("list", length(core.fams))
  for (j in seq_along(core.fams)) {

    fam <- core.fams[[j]]
    curr.fam.name <- names(core.fams)[j]

    # Identify all possible genotype combinations for the core families
    founders <- which(fam$mother == 0 & fam$father == 0)
    genos <- matrix(0, length(founders), nrow(fam))
    for (i in seq_along(founders)) {
      genos[i, founders[i]] <- 1
    }

    expand <- function(genos) {
      # Add more possible genotype combinations by adding children of carriers
      # to the list of carriers
      f <- function(g) {
        carriers <- which(g == 1)
        children <- which(fam$mother %in% carriers | fam$father %in% carriers)
        children <- setdiff(children, carriers)
        if (length(children) == 0) {
          extra <- matrix(0, 0, nrow(fam))
        } else {
          extra <- matrix(g, length(children), nrow(fam), byrow = TRUE)
          for (i in seq_along(children)) {
            extra[i, children[i]] <- 1
          }
        }
        extra
      }

      extra.genos <- vector("list", nrow(genos))
      for (i in 1:nrow(genos)) {
        extra.genos[[i]] <- f(genos[i, ])
      }

      # Collate into a matrix
      genos <- do.call(rbind, extra.genos)

      # Remove duplicates
      binary <- function(g) {
        sum(g * 2^((nrow(fam) - 1):0))
      }
      unique.number <- apply(genos, 1, binary)
      keep <- !duplicated(unique.number)
      genos <- genos[keep, ]
    }

    genos.withgivencarriers <- genos
    n.old.genos <- 0
    num.carriers <- 1

    while (n.old.genos < nrow(genos)) {
      num.carriers <- num.carriers + 1
      n.old.genos <- nrow(genos)
      genos.withgivencarriers <- expand(genos.withgivencarriers)
      genos <- rbind(genos, genos.withgivencarriers)
    }
    dim(genos)

    # Create matrix of transmission probablities
    trans.matrix <- create_trans_matrix()

    trans.prob <- function(g) {
      gm <- g[fam$mother]
      gd <- g[fam$father]
      go <- g[fam$mother != 0 & fam$father != 0]
      y <- 1 + 9 * go + 3 * gm + gd
      p <- trans.matrix$p[y]
      prod(p)
    }

    afreq <- 0.01
    p <- apply(genos, 1, trans.prob)
    c(sum(p), length(founders)) # Should be the same
    genos <- rbind(rep(0, nrow(fam)), genos)
    p <- c((1 - afreq)^2, 2 * afreq * (1 - afreq) * p)
    p <- p / sum(p)
    genos <- data.frame(genos, p)
    names(genos)[1:nrow(fam)] <- as.character(fam$ID.orig)

    results[[j]] <- genos

  }
  names(results) <- names(core.fams)
  results
}


