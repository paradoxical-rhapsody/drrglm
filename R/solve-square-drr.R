#' @title Doubly Regularized Square Loss
#' @description
#' Solve the following problem
#' \deqn{
#'     \min_{L, S} \Big\lbrace \frac{1}{2} \| Z - (L+S) \|_F^2 +
#'      \lambda_1 \|L\|_* + \lambda_2 \| S \|_1 \Big\rbrace.
#' }
#' 
#' @param Z Numeric matrix.
#' @param lambda1 \eqn{\lambda_1}.
#' @param lambda2 \eqn{\lambda_2}.
#' @param S.diag.penalize Whether to penalize the diagonal elements of `S`.
#'  `S.diag.penalize=FALSE` is used in factor model
#'  (See [tune_drr_factor_model], default `TRUE`).
#' @param tol Convergence tolerance for \eqn{L}.
#' @param maxIter Maximal step of iterations for \eqn{L}.
#' 
#' @return `list(L, S, Lsv, iter, isConvergent, lambda1, lambda2)`.
#' 
#' @examples
#' set.seed(2023)
#' r0 <- 9
#' c0 <- 7
#' Z <- matrix(runif(r0*c0), r0)
#' lambda1 <- 0.5
#' lambda2 <- 0.1
#' S.diag.penalize <- TRUE
#' 
#' tol <- 1e-3
#' maxIter <- 200
#' 
#' system.time( egg <- solve_square_drr(Z, lambda1, lambda2, S.diag.penalize, tol, maxIter) )
#' 
#' @noRd
solve_square_drr <- function(Z, lambda1, lambda2, S.diag.penalize, tol, maxIter){
    stopifnot( length(lambda1) == 1 && length(lambda2) == 1 )
    stopifnot( lambda1 > 0 && lambda2 > 0 )

    Ltmp <- cxx_solve_huber_nuclear(Z, lambda2, lambda1, tol, maxIter)
    S <- solve_square_lasso(Z - Ltmp[["L"]], lambda2, S.diag.penalize)

    list(
        L = Ltmp[["L"]],
        S = S,
        Lsv = Ltmp[["svals"]],
        iter = Ltmp[["iter"]],
        isConvergent = Ltmp[["isConvergent"]],
        lambda1 = lambda1,
        lambda2 = lambda2
    )
}
