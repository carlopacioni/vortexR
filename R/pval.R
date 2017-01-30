#' Calculates p-values from z-values
#'
#' \code{pval} calculates one-tailed p values from a vector that contains
#'   z-values and it is geenrally used internally.
#'
#' @param x z-values
#' @return A numeric vector of length equal to length(x)
#' @export
#' @examples
#'  z <- c(1.645, 1.96, 3.09)
#'  pval(z)
pval <- function(x) pnorm(abs(x), lower.tail=FALSE)
