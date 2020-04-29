# Change the geno probs in typed.genos to reflect a given allele freq
adjust_geno_probs <- function(typed.genos, afreq = 0.01) {
  # Also, condition the likelihood on the event that all of the founders
  # together carry at most one mutation (otherwise the likelihood for the null
  # model won't match the likelihood calculated below). So afreq needs to be
  # small for this conditional likelihood to be approximately the same as the
  # likelihood function for the full segregation analysis model.
  for (i in seq_along(typed.genos)) {
    tg <- typed.genos[[i]]
    nfounders <- sum(tg$p) / 2
    if (abs(nfounders - round(nfounders)) > 1e-5) {
      warning("Non-integer number of founders for family ", i, call. = FALSE)
    }
    nfounders <- round(nfounders)
    tg$p <- tg$p * afreq * (1 - afreq)^(2 * nfounders - 1)
    f <- function(x) all(x[-1] == 0)
    j <- which(apply(tg, 1, f))
    prob <- (1 - afreq)^(2 * nfounders)
    if (length(j) == 0) {
      x <- c(prob, rep(0, ncol(tg) - 1))
      tg <- rbind(x, tg)
    }
    if (length(j) == 1) {
      tg$p[j] <- tg$p[j] + prob
    }
    if (length(j) > 1) {
      j <- j[1]
      tg$p[j] <- tg$p[j] + prob
      warning("Repeated genotype combinations for family ", i, call. = FALSE)
    }
    # Condition on the event that the genotype combinations in tg are approx
    # exhaustive
    tg$p <- tg$p / sum(tg$p)
    typed.genos[[i]] <- tg
  }
  typed.genos
}

