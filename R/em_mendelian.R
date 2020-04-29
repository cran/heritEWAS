#' @importFrom stats dnorm sd
em_mendelian <- function(typed.genos, M, ncores = 1) {

  indices.typed.genos <- find_geno_index(typed.genos, M)

  update.theta <- function(theta, mval) {
    # Update of the parameter vector theta=c(mu0,sd0,mu1,sd1) with the EM
    # algorithm
    # Here our likelihood is a Gaussian mixture, with different Gaussian
    # densities for carriers (params mu0 and sd0) and non-carriers (params mu1
    # and sd1)
    # of a rare mutation which affects methylation values at a probe with
    # M-values mval.
    # This is the standard algorithm but modified to take non-independence of y
    # into account, where:
    #   * the data is x = (x1, ..., xn) is the vector of M-values mval
    #   * the hidden data y = (y1, ..., yn) gives the true genotypes of each
    #     person

    # Decode the parameter vector theta
    mu0 <- theta[1]
    sd0 <- theta[2]
    mu1 <- theta[3]
    sd1 <- theta[4]

    # Calculate q where
    # q[i] = q_{i1}^t in the notation of my notes
    #      = probability that person i (matching mval) is a carrier based on
    #        the M-values and previous values of the parameters
    y <- c(dnorm(mval, mu0, sd0), dnorm(mval, mu1, sd1))
    p01 <- matrix(y, length(mval), 2)
    f <- function(j) {
      # j is the family number
      datg <- typed.genos[[j]]
      p01.fam <- p01[indices.typed.genos[[j]], ]
      if (ncol(datg) == 2) p01.fam <- matrix(p01.fam, 1, )
      p <- datg$p
      for (k in 2:ncol(datg)) {
        indices <- datg[, k] + 1
        p <- p * p01.fam[k - 1, indices]
      }
      p <- p / sum(p)
      g <- function(k) {
        carriers <- (datg[, k] == 1)
        sum(p[carriers])
      }
      q.famj <- sapply(2:ncol(datg), g)
      q.famj
    }
    q.fams <- lapply(1:length(typed.genos), f)
    q <- numeric(length(mval))
    for (j in 1:length(typed.genos)) {
      q[indices.typed.genos[[j]]] <- q.fams[[j]]
    }

    # Calculate the updated estimates
    w0 <- (1 - q) / sum(1 - q)
    w1 <- q / sum(q)
    mu0 <- sum(w0 * mval)
    sd0 <- sqrt(sum(w0 * (mval - mu0)^2))
    mu1 <- sum(w1 * mval)
    sd1 <- sqrt(sum(w1 * (mval - mu1)^2))

    # Output
    theta <- c(mu0, sd0, mu1, sd1)
    theta
  }

  calc.ll.em <- function(theta, mval) {
    # Calculate the likelihood whcih is maximized by the EM algorithm
    # for given parameters theta=c(mu0,sd0,mu1,sd1) and M-values at a
    # particular probe mval

    # Decode the parameter vector theta
    mu0 <- theta[1]
    sd0 <- theta[2]
    mu1 <- theta[3]
    sd1 <- theta[4]

    # Calculate P(x | theta)
    y <- c(dnorm(mval, mu0, sd0), dnorm(mval, mu1, sd1))
    p01 <- matrix(y, length(mval), 2)
    f <- function(j) {
      # j is the family number, this calculates the log-likelihood for family j
      datg <- typed.genos[[j]]
      p01.fam <- p01[indices.typed.genos[[j]], ]
      if (ncol(datg) == 2) p01.fam <- matrix(p01.fam, 1, )
      p <- datg$p
      for (k in 2:ncol(datg)) {
        indices <- datg[, k] + 1
        p <- p * p01.fam[k - 1, indices]
      }
      ll.famj <- log(sum(p))
      ll.famj
    }
    p.fams <- sapply(1:length(typed.genos), f)

    # Output
    ll <- sum(p.fams)
    ll
  }

  em.maxlik <- function(i, tol = 1e-3, max.count = 200, theta.start = NULL) {
    # Use the EM algorithm to produce ML estimates for probe i

    # Initialize
    mval <- as.numeric(M[i, ])
    n <- length(mval)
    mu0 <- mean(mval)
    sd0 <- sd(mval) * sqrt((n - 1) / n)
    mu1 <- mu0
    sd1 <- sd0
    theta.null <- c(mu0, sd0, mu1, sd1)
    if (is.null(theta.start)) {
      theta.start <- theta.null
    }
    theta <- theta.start

    # Iterate
    continue <- TRUE
    count <- 0
    while (continue) {
      count <- count + 1
      theta.new <- update.theta(theta, mval)
      continue <- all(!is.na(theta.new)) & any(abs(theta.new - theta) > tol) &
        (count < max.count)
      theta <- theta.new
    }

    # Calculate likelihoods
    # This should be called ll.null, but this name kept for compatibility with
    # other scripts
    ll.start <- calc.ll.em(theta.null, mval)
    ll.end <- calc.ll.em(theta, mval)

    # Output
    c(ll.start, ll.end, theta.null[1:2], theta, count)
  }

  # Choose three stating values and choose the best EM output
  em.maxlik.best <- function(i) {

    # Three sets of starting values for EM
    theta_start_values <- list(
      NULL, c(-2, 1, 2, 1), c(2, 1, -2, 1)
    )

    out <- vector("list", length = 3)
    for (j in seq_along(out)) {
      out[[j]] <- em.maxlik(i, theta.start = theta_start_values[[j]])
    }
    out <- matrix(unlist(out, use.names = FALSE), nrow = 3, byrow = TRUE)

    #max_out <- out[which.max(out[, 2]), ]

    out_noinf <- rbind(out[out[, 2] != Inf, ])
    max_out <- out_noinf[which.max(out_noinf[, 2]), ]

    #if (length(max_out) == 0 | max_out[2] == Inf) {
    if (length(max_out) == 0) {
      max_out <- out[1, ]
      max_out[-c(1, 3, 4)] <- NA
    }
    max_out
  }

  nprobe <- nrow(M)
  if (ncores == 1) {
    out <- sapply(1:nprobe, em.maxlik.best)
  } else {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    out <- parallel::parSapply(cl, 1:nprobe, em.maxlik.best)
  }
  out <- data.frame(probe = rownames(M)[1:nprobe], as.data.frame(t(out)))
  names(out) <- c("probe", "ll.null", "ll.mendel", "mu.null", "sd.null",
                  "mu0", "sd0", "mu1", "sd1", "count")
  out

}

