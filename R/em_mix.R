#' @importFrom stats dnorm sd
em_mix <- function(typed.genos, M, ncores = 1) {

  update.theta <- function(theta, mval, afreq = 0.01) {
    # Update of the parameter vector theta=c(mu0,sd0,mu1,sd1) with the standard
    # EM algorithm for Gaussian mixtures (with two groups, 0 and 1) but with
    # the prior group membership probabilities fixed.
    # Group membership doesn't have any genetic meaning and is assumed to be
    # independent across individuals and probes.
    # We use the standard EM algorithm where:
    #   * the data is x = (x1, ..., xn) is the vector of M-values mval at a
    #     given probe
    #   * the hidden data y = (y1, ..., yn) gives the group membership of
    #     each person

    # Decode the parameter vector theta
    mu0 <- theta[1]
    sd0 <- theta[2]
    mu1 <- theta[3]
    sd1 <- theta[4]

    # Group membership probabilities based on x and theta
    # (i.e. theta_t, from the previous step)
    # These should be estimated (and so included in theta) if the carrier
    # probbaility was estimated
    # for the Mendelian model.  Otherwise they should be fixed to the same
    # values of the group (i.e. carrier) probabilities used in the Mendelian
    # model.
    alpha0 <- 1 - afreq * (1 - afreq)
    alpha1 <- afreq * (1 - afreq)

    # Calculate q = P(group=1 | x_i, theta_t)
    p0 <- alpha0 * dnorm(mval, mu0, sd0)
    p1 <- alpha1 * dnorm(mval, mu1, sd1)
    q <- p1 / (p0 + p1)
    q[p1 == Inf & p0 != Inf] <- 1

    # Calculate the updated estimates
    w0 <- (1 - q) / sum(1 - q)
    w1 <- q / sum(q)
    mu0 <- sum(w0 * mval)
    sd0 <- sqrt(sum(w0 * (mval - mu0)^2))
    mu1 <- sum(w1 * mval)
    sd1 <- sqrt(sum(w1 * (mval - mu1)^2))
    c(mu0, sd0, mu1, sd1)

    # Output
    theta <- c(mu0, sd0, mu1, sd1)
    theta
  }

  calc.ll.em <- function(theta, mval, afreq = 0.01) {
    # Calculate the likelihood whcih is maximized by the EM algorithm
    # for given parameters theta=c(mu0,sd0,mu1,sd1) and M-values at a
    # particular probe mval

    # Decode the parameter vector theta
    mu0 <- theta[1]
    sd0 <- theta[2]
    mu1 <- theta[3]
    sd1 <- theta[4]

    # Group membership probabilities based on x and theta
    # (i.e. theta_t, from the previous step)
    # These should be estimated (and so included in theta) if the carrier
    # probbaility was estimated for the Mendelian model.
    # Otherwise they should be fixed to the same values of the group
    # (i.e. carrier) probabilities used in the Mendelian model.
    alpha0 <- 1 - afreq * (1 - afreq)
    alpha1 <- afreq * (1 - afreq)

    # Calculate P(x | theta)
    # where P(xi | theta) = sum_yi P(xi | yi, theta)*P(yi | theta)
    #                     = sum_yi f_yi(xi)*alpha_yi
    f0 <- dnorm(mval, mu0, sd0)
    f1 <- dnorm(mval, mu1, sd1)
    p <- alpha0 * f0 + alpha1 * f1

    # Output
    ll <- sum(log(p))
    ll
  }

  em.maxlik <- function(i, tol = 1e-3, max.count = 200, theta.start = NULL) {
    # Use the EM algorithm to produce ML estimates for probe i

    # Initialize
    mval <- as.numeric(M[i, ])
    n <- length(mval)
    mu0 <- mean(mval)
    sd0 <- sd(mval) * sqrt((n - 1) / n)
    theta.null <- c(mu0, sd0, mu0, sd0)
    if (is.null(theta.start)) {
      # If the densities are the same for the two groups then theta never
      # changes, I think
      mu1 <- mu0 + 0.01
      sd1 <- sd0 * 1.01
      theta.start <- c(mu0, sd0, mu1, sd1)
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
    ll.null <- calc.ll.em(theta.null, mval)
    ll.end <- calc.ll.em(theta, mval)

    # Output
    c(ll.null, ll.end, theta.null[1:2], theta, count)
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
  names(out) <- c("probe", "ll.null", "ll.mix", "mu.null", "sd.null",
                  "mu0", "sd0", "mu1", "sd1", "count")
  out

}
