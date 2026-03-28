#' @title Simulation: Setups in Regularized Matrix Regression
#' 
#' @param N Sample size.
#' @param p1 Row dimension.
#' @param p2 Column dimension.
#' @param r0 Rank.
#' @param s0 signal proportion.
#' @param family `gaussian` or `binomial`.
#' 
#' @return `list( C0, train=list(x, y), test=list(x, y) )`.
#' 
#' @references
#' 
#' Zhou Hua and Li Lexin (2014). Regularized Matrix Regression. Journal of the Royal Statistical Society (Series B). <doi:10.1111/rssb.12031>
#' 
#' @examples
#' N <- 100
#' p1 <- 64
#' p2 <- 64
#' r0 <- 1
#' s0 <- 0.05
#' family <- 'gaussian'
#' 
#' set.seed(2026)
#' DT <- simu_zhouandli2014(N, p1, p2, r0, s0, family)
#' 
#' @export
simu_zhouandli2014 <- function(N, p1, p2, r0, s0, family){
    stopifnot( s0 > 0 && s0 < 1)
    stopifnot( family %in% c("gaussian", "binomial") )

    B <- tcrossprod(
        matrix(rbinom(p1*r0, 1, sqrt(1 - (1-s0)^(1/r0))), p1),
        matrix(rbinom(p2*r0, 1, sqrt(1 - (1-s0)^(1/r0))), p2)
    )
    stopifnot( p1 == NROW(B) )
    stopifnot( p2 == NCOL(B) )

    set_x_and_y <- function(family, N){
        x <- array(rnorm(p1*p2*N), c(p1, p2, N))

        eta <- apply(x, 3, function(xi) sum(B * xi))
        if (family == "gaussian") {
            y <- eta + rnorm(N, sd=1)
        }
        if (family == "binomial") {
            y <- rbinom(N, 1, binomial()$linkinv(eta))
        }

        stopifnot( length(y) == dim(x)[3] )
        return( list(x=x, y=y) )

    }

    list(
        C0 = B,
        train = set_x_and_y(family, N),
        test  = set_x_and_y(family, 1000)
    )
}
