check_dat <- function(dat) {

  # Check if the input is a data.frame
  if (class(dat) != "data.frame") {
    #stop(paste0(deparse(substitute(dat)), " is not a data.frame"))
    stop(sprintf("%s is not a data.frame", deparse(substitute(dat))),
         call. = FALSE)
  }

  # Check if all required columns are here
  column_names <- c("family", "indiv", "mother", "father", "typed")
  w <- which(!column_names %in% names(dat))
  if (!all(column_names %in% names(dat))) {
    stop(paste("Cannot find column", paste(sQuote(column_names[w]),
                                            collapse = ", ")), call. = FALSE)
  }

  # Get columns
  family <- dat$family
  indiv <- dat$indiv
  mother <- dat$mother
  father <- dat$father
  typed <- dat$typed

  # Check if any missing values in column 'family', 'indiv' and 'typed'
  na_cols <- sapply(list(family, indiv, typed), anyNA)
  if (sum(na_cols) >= 1) {
    cols <- sQuote(column_names[c(1, 2, 5)][na_cols])
    stop(paste("Column", paste(cols, collapse = ", "),
               "should not have missing values"), call. = FALSE)
  }

  # Check if only 0 or 1 are present in column 'typed'
  only_zero_one <- isTRUE(all.equal(sort(unique(typed)), 0:1))
  if (!only_zero_one) {
    stop("Column 'typed' should only contain values 0 and 1", call. = FALSE)
  }

  # Check if 'indiv' column has unique values
  if (sum(duplicated(indiv)) >= 1) {
    stop("Column 'indiv' should be unique for each person", call. = FALSE)
  }

}

