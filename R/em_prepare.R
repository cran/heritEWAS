# Delete any columns of M that aren't included in the pedigree dataset
filter_M <- function(typed.genos, M) {
  typedIDs <- unlist(lapply(typed.genos, names))
  typedIDs <- setdiff(typedIDs, "p")
  keep <- colnames(M) %in% typedIDs
  M <- M[, keep]
  M
}

# Find the columns of M corresponding to the columns of each element of
# typed.genos
find_geno_index <- function(typed.genos, M) {
  f <- function(l) {
    x <- colnames(l)[-1]
    index <- numeric(length(x))
    for (i in seq_along(x)) {
      index[i] <- which(colnames(M) == x[i])[1]
    }
    index
  }
  lapply(typed.genos, f)
}
