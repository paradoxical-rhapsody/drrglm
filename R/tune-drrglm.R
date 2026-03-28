#' @name tune.drrglm
#' @title Doubly Regularized Regression for Matrix-Variate Data
#' @description Tune the \eqn{(\lambda_1, \lambda_2)} in [drr].
#' 
#' @param iniParas Initial values returned by [ini_paras]. See examples below.
#' @param x Numeric array.
#' @param y Numeric vector.
#' @param lambda1.factor Factor on initial guess\eqn{\lambda_1^2}
#'  (By default `1.3^seq(-3, 3, 0.2)`).
#' @param maxCard The maximal cardinality of \eqn{S}.
#' @param tol Convergence tolerance.
#' @param maxIter Maximal step of iterations.
#' 
#' @return List.
#' 
#' @examples
#' p1 <- 10
#' p2 <- 9
#' rank0 <- 3
#' S.type <- "rand"
#' coefs <- simu_reg_coefs(p1, p2, rank0, S.type)
#' 
#' N <- 500
#' family <- "gaussian"
#' L0 <- coefs[["L"]]
#' S0 <- coefs[["S"]]
#' C0 <- L0 + S0
#' DT <- simu_reg_data(family, N, C0)
#' 
#' x <- DT[["train"]][["x"]]
#' y <- DT[["train"]][["y"]]
#' 
#' lambda1.factor <- 1.1^(0:1)
#' maxRank <- 5
#' 
#' system.time( iniParas <- ini_paras(x, y, family, maxRank) )
#' system.time( drrObj <- tune_drrglm(iniParas, x, y, lambda1.factor) )
#' 
#' @export
tune_drrglm <- function(iniParas, x, y, lambda1.factor=1.3^seq(-3, 3, 0.2), 
                        maxCard=30, tol=1e-3, maxIter=500){
    stopifnot( length(dim(x)) == 3 )
    stopifnot( length(y) == dim(x)[3] )

    valuetol <- min(0.1*tol, 1e-3)

    rDim <- NROW(x)
    cDim <- NCOL(x)
    N <- dim(x)[3]

    family <- get("family", iniParas)
    grad <- get("grad", iniParas)
    lipschitz <- get("lipschitz", iniParas)

    sigma0 <- get("sigma0", iniParas)
    rankStar <- get("rankStar", iniParas)
    rankStarHat <- get("rankStarHat", iniParas)
    lambda.star <- get("lambda.star", iniParas)
    L0.star <- get("L0.star", iniParas)

    lambda1Set <- sort(lambda1.factor, decreasing=TRUE) * lambda.star
    stopifnot( length(lambda1Set) == length(lambda1.factor))

    loss0 <- Inf
    lik.bic0 <- Inf
    loss.star <- Inf
    C0 <- L0 <- matrix(0.0, NROW(x), NCOL(x))
    for (lambda1 in lambda1Set) {
        nrrobj <- nrrglm(x, y, family, lambda1, L0, tol, maxIter=500, FALSE)
        L0 <- nrrobj[["L"]]
        eta.L0 <- apply(x, 3, function(xn) sum(xn * L0))
        stopifnot( dim(L0) == c(NROW(x), NCOL(x)) )

        lambda2Set <- glmnet(
            x = matrix(x, N, byrow=TRUE),
            y = y,
            family = family,
            offset = eta.L0,
            standardize = FALSE,
            intercept = FALSE
        )[["lambda"]]

        lambda2Set <- lambda2Set[seq( min(length(lambda2Set), maxCard) )]

        for (lambda2 in lambda2Set) {
            drrobj <- drrglm(x, y, family, lambda1, lambda2, C0, tol, maxIter, verbose=FALSE)
            stopifnot( dim(drrobj$L) == c(NROW(x), NCOL(x)) )
            stopifnot( dim(drrobj$S) == c(NROW(x), NCOL(x)) )

            C0 <- drrobj[["L"]] + drrobj[["S"]]
            rankL <- sum(abs(drrobj[["Lsv"]]) >= valuetol)
            cardS <- sum(abs(drrobj[["S"]]) >= valuetol)

            dof <- cardS
            eta.train <- drop( apply(x, 3, function(xn) sum(xn*C0)) )
            glmObj <- glm(y ~ 0, family=family, offset=eta.train)

            lik.bic <- -2*logLik(glmObj) + dof * log(N)
            if (lik.bic < lik.bic0) { 
                drr.lik.bic <- drrobj
                lik.bic0 <- lik.bic
            }

        }
    }

    class(drr.lik.bic) <- "drr"
    return(drr.lik.bic)
}
