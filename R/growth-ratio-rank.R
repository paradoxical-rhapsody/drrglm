#' @title Identify Rank using Growth-Ratio Criterion
#' @description
#' Let \eqn{V(k) = \sum_{j=k+1}^r \sigma_j (\hat{L})},
#' \deqn{GR(k) = \frac{
#'      \log(1 + \sigma_k(\hat{L}) / V(k))}{
#'      \log(1 + \sigma_{k+1}(\hat{L})/V(k+1))}
#' }
#' 
#' @param sv Vector of singular values.
#' @param tol The singular values smaller than tol is treated as zero.
#' 
#' @return An integer of rank.
#' 
#' @examples
#' set.seed(2025)
#' A1 <- matrix(rnorm(10*3), 10)
#' A2 <- matrix(rnorm(12*3), 12)
#' sv <- svd(A1 %*% t(A2))$d

#' ( egg <- get_gr_rank(sv) )
#' 
#' @noRd
get_gr_rank <- function(sv, tol=1e-4){
    stopifnot( all(sv >= 0) )

    if ( all(sv < tol) ) 
        return(0)

    n <- length(sv)

    tailsum <- function(k)
        log(1 + sv[k]/sum(sv[(k+1):n])) / log(1+sv[k+1]/sum(sv[(k+2):n]))

    GR <- sapply(1:(n-2), tailsum)

    return( which.max(GR) )
}
