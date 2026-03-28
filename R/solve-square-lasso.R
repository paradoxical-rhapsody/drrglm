#' @title Soft Thresholding
#' @description
#' \deqn{
#'  \begin{align*}
#'     \hat{S}
#'       &= \arg \min_{S} \Big\lbrace \frac{1}{2} \| Z - S \|_F^2 + \lambda \| S \|_1 \Big\rbrace \\\\\\
#'       &= \Big[ \text{sign}(Z_{ij}) \max\big( |Z_{ij}| - \lambda, 0 \big) \Big].
#'  \end{align*}
#' }
#' 
#' @param Z Numeric matrix.
#' @param lambda Tuning parameter \eqn{\lambda}.
#' @param S.diag.penalize Whether to penalize the diagonal elements of `S`.
#'  `S.diag.penalize=FALSE` is used in factor model
#'  (See [tune_drr_factor_model], default `TRUE`).

#' 
#' @return Matrix `S`.
#' 
#' @examples
#' set.seed(2023)
#' r0 <- 7
#' c0 <- 5
#' Z <- matrix(rnorm(r0*c0), r0)
#' lambda <- 1.5
#' S.diag.penalize <- FALSE
#' 
#' system.time( S <- solve_square_lasso(Z, lambda, S.diag.penalize) )
#' 
#' @noRd
solve_square_lasso <- function(Z, lambda, S.diag.penalize){
    stopifnot( length(lambda) == 1 )
    stopifnot( lambda > 0 )
    
    egg <- sign(Z) * pmax(abs(Z) - lambda, 0.0)
    if (isFALSE(S.diag.penalize))
        diag(egg) <- diag(Z)

    return(egg)
    # objVal <- 0.5 * norm(Z-S, "F")^2 + lambda*sum(abs(S))
    # return(list(S=S, objVal=objVal))
}
