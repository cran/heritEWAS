#' Simulated data on 20 families.  
#' 
#' A dataset giving the relationship structure of 20 families 
#' and phenotypic data on the family members
#'
#' @format A data frame with 288 rows (corresponding to persons) and 8 variables:
#' \describe{
#'   \item{family}{an identifier for the person's family}
#'   \item{indiv}{an identifier (ID) for the person}
#'   \item{mother}{the individual ID of the person's mother}
#'   \item{father}{the individual ID of the person's father}
#'   \item{sex}{the person's sex (M = male, F = female)}
#'   \item{aff}{the person's affected status (1 = case, 0 = control)}
#'   \item{age}{the person's age, in years}
#'   \item{typed}{a flag indicating if methylation data is available (1 = available, 0 = unavailable)}
#' }
#' @source Simulated
"ped"