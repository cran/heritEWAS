trim_pedigree <- function(fam, reorder = TRUE) {

  # Convert to fam IDs to 1,...,n
  fam <- convert_IDs(cbind(family = 1, fam), convert.IDs.numeric = TRUE)[, -1]

  still.updating <- TRUE
  while (still.updating) {
    still.updating <- FALSE

    # First, restrict to typed people and their ancestors
    ID.list <- fam$indiv[fam$typed == 1]
    keep <- fam$indiv %in% ancestors(ID.list, fam)
    if (any(!keep)) still.updating <- TRUE
    fam <- fam[keep, ]

    # Second, remove pairs of founders who link to typed persons always through
    # one of their children
    ID.list <- fam$indiv[fam$typed == 1]
    ap <- ancestor.paths(ID.list, fam)
    founders <- fam$indiv[fam$mother == 0 & fam$father == 0]
    while (length(founders) > 0) {
      f <- founders[1]
      child <- overlap(f, ap)
      if (child == -1) {
        founders <- founders[-1]
      } else {
        parents <- c(fam$mother[fam$indiv == child],
                     fam$father[fam$indiv == child])
        if (all(parents %in% founders)) {
          still.updating <- TRUE
          fam$mother[fam$mother %in% parents] <- 0
          fam$father[fam$father %in% parents] <- 0
          fam <- fam[!(fam$indiv %in% parents), ]
        }
        founders <- setdiff(founders, parents)
      }
    }
  }

  if (reorder) {
    fam <- convert_IDs(cbind(family = 1, fam), convert.IDs.numeric = TRUE)[, -1]
  }
  fam
}

