# Return ID.list and all of their ancestors in fam
ancestors <- function(ID.list, fam) {
  x <- c()
  y <- ID.list
  while (length(x) < length(y)) {
    x <- y
    flag <- fam$indiv %in% x
    y <- unique(c(x, fam$mother[flag], fam$father[flag]))
    y <- setdiff(y, 0)
  }
  y
}

ancestors2 <- function(ID.list, fam) {
  x <- c()
  y <- ID.list
  while (length(x) < length(y)) {
    x <- y
    flag <- fam$indiv %in% x
    y <- unique(c(x, fam$mother[flag], fam$father[flag]))
    y <- setdiff(y, "")
    y <- setdiff(y, 0)
  }
  y
}

# Return ID.list and all of their descendents in fam
descendents <- function(ID.list, fam) {
  x <- c()
  y <- ID.list
  while (length(x) < length(y)) {
    x <- y
    flag <- (fam$mother %in% x) | (fam$father %in% x)
    y <- unique(c(x, fam$indiv[flag]))
  }
  y
}

# Return a list of all ancenstor paths in fam starting at any person in
# ID.list, where an ancestor path is a vector p of IDs so that p[i+1] is a
# parent of p[i] for all i
ancestor.paths <- function(ID.list, fam) {
  x <- c()
  y <- ID.list
  while (length(x) < length(y)) {
    x <- y
    y <- c()
    for (i in 1:length(x)) {
      p <- x[[i]]
      last <- fam$indiv == p[length(p)]
      if (fam$mother[last] == 0 & fam$father[last] == 0) y <- c(y, list(p))
      if (fam$mother[last] != 0) y <- c(y, list(c(p, fam$mother[last])))
      if (fam$father[last] != 0) y <- c(y, list(c(p, fam$father[last])))
    }
  }
  y
}

# If all ancestor paths p ending at founder f pass through a child of f then
# return the child's ID, otherwise return -1
overlap <- function(founder, ap) {
  g <- function(p) {
    return(p[length(p)] == founder)
  }
  keep <- sapply(ap, g)
  apf <- ap[keep]
  h <- function(p) {
    if (length(p) == 1) {
      return(-1)
    }
    return(p[length(p) - 1])
  }
  children <- unique(unlist(lapply(apf, h)))
  if (any(children == -1)) {
    return(-1)
  }
  if (length(children) > 1) {
    return(-1)
  }
  return(children[1])
}


# Generate combinations of possible genotypes
geno.comb <- function(factors, m) {
  X <- data.frame(v = factors)
  if (m == 1) {
    return(X)
  }
  for (i in 2:m) {
    # i <- 2
    Y <- c()
    for (f in factors) {
      Xf <- data.frame(v = rep(f, nrow(X)), X)
      Y <- rbind(Y, Xf)
    }
    X <- Y
  }
  names(X) <- paste0("V", 1:ncol(X))
  X
}

# Generate possible genotypes for typed people, assume only one founder carrier
generate_typed_geno <- function(fam) {

  # Find the founders, typed people and their intersection
  founders <- fam$indiv[fam$mother == 0 & fam$father == 0]
  typed <- fam$indiv[fam$typed == 1]
  common <- intersect(founders, typed)

  # Check if all typed people have a common (founder) ancestor
  common.ancestor <- FALSE
  for (i in 1:length(founders)) {
    if (all(typed %in% descendents(founders[i], fam))) {
      common.ancestor <- TRUE
    }
  }

  if (common.ancestor) {
    typed.geno <- geno.comb(0:1, length(typed))
  } else {
    typed.geno <- c()
    for (i in seq_along(founders)) {
      td <- intersect(typed, descendents(founders[i], fam))
      if (length(td) == 1) {
        extra <- data.frame(matrix(0, 2, length(typed)))
        extra[, typed == td] <- c(0, 1)
        typed.geno <- rbind(typed.geno, extra)
      } else if (length(td) > 1) {
        geno.td <- geno.comb(0:1, length(td))
        extra <- data.frame(matrix(0, nrow(geno.td), length(typed)))
        extra[, typed %in% td] <- geno.td
        typed.geno <- rbind(typed.geno, extra)
      }
    }
    # Remove duplicate rows
    y <- numeric(nrow(typed.geno))
    for (i in 1:ncol(typed.geno)) {
      y <- y + typed.geno[, i] * 2^(i - 1)
    }
    typed.geno <- typed.geno[!duplicated(y), ]
  }
  typed.geno
}

# Create matrix of transmission probablities
create_trans_matrix <- function() {

  gamete <- function(off, par) {
    ifelse(par == 1, 0.5, as.numeric(par == 2 * off))
  }

  trans.matrix <- c()
  for (go in 0:2) {
    for (gm in 0:2) {
      for (gd in 0:2) {
        p <- -1
        if (go == 0) p <- gamete(0, gm) * gamete(0, gd)
        if (go == 2) p <- gamete(1, gm) * gamete(1, gd)
        if (go == 1) {
          p <- gamete(0, gm) * gamete(1, gd) +
            gamete(1, gm) * gamete(0, gd)
        }
        trans.matrix <- rbind(trans.matrix, c(go, gm, gd, p))
      }
    }
  }
  trans.matrix <- data.frame(trans.matrix)
  names(trans.matrix) <- c("go", "gm", "gd", "p")
  trans.matrix
}

