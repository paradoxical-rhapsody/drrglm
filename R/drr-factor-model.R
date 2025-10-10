#' @title Simulating Data for Factor Model with Correlated Noise
#' @description DRR for factor model.
#' 
#' @param y Numeric matrix of \eqn{N \times p}.
#' @param lambda1.factor Factor on \eqn{\lambda_1^8}
#'  (By default `1.1^(-2:2)`).
#' @param S.diag.penalize Whether to penalize the diagonal elements of `S` (default `FALSE`).
#' @param tol Convergence tolerance.
#' @param maxIter Maximal step of iterations.
#' 
#' @return List.
#' 
#' @examples
#' set.seed(2025)
#' N <- 500
#' p <- 30
#' r <- 10
#' noise <- "diag"
#' paras <- simu_factor_model_paras(p, r, noise)
#' y <- simu_factor_model_data(N, paras[["B"]], paras[["S"]])
#' lambda1.factor <- 1.1^(-1:1)
#' system.time( egg <- tune_drr_factor_model(y, lambda1.factor) )
#' 
#' @export
tune_drr_factor_model <- function(y, lambda1.factor=1.1^(-2:2), S.diag.penalize=FALSE, 
                tol=1e-3, maxIter=500){
    valuetol <- 1e-6
    N <- NROW(y)
    p <- NCOL(y)

    Z <- cov(y)
    Z.norm <- sum( (Z[upper.tri(Z, diag=TRUE)])^2 )
    Zsvd <- svd(Z)
    d0 <- Zsvd[["d"]]

    ER <- d0[-length(d0)] / d0[-1]
    rstar <- which.max(ER) # min( which.max(ER), 20 )
    lambda1.star <- 0.01*d0[rstar] + 0.99*d0[rstar + 1]
    stopifnot( length(rstar) == 1 )
    stopifnot( length(lambda1.star) == 1 )

    ## NRR estimator
    nrrd0 <- pmax(d0[1:rstar] - lambda1.star, 0.0)
    nrrL <- Zsvd[["u"]][, 1:rstar, drop=FALSE] %*% 
        diag(nrrd0, nrow=length(nrrd0)) %*%
        t(Zsvd[["u"]][, 1:rstar, drop=FALSE])
    stopifnot( nrrd0 > 0 )
    stopifnot( length(nrrd0) == rstar )

    ## DRR estimator
    maxCard <- min(p*(p-1)/2, 60)
    stopifnot( maxCard <= p*(p-1)/2 )   # max number of nonzero off-diagonal entries
    lambda1Set <- lambda1.factor * lambda1.star
    
    cnames <- c("lambda1", "lambda2", "rankL", "cardS", "loss", "ebic", "isConvergent")
    performance <- matrix(NA_real_, length(lambda1Set)*maxCard, length(cnames))
    colnames(performance) <- cnames

    k <- 1
    L0 <- nrrL
    loss0 <- ebic0 <- Inf
    ebic.gamma <- max(0.0, 1.0 - log(N)/(2.0*log(p*(p+1)/2)))
    for (lambda1 in lambda1Set) {
        Ltmp <- Zsvd[["u"]] %*% 
            diag(pmax(d0-lambda1, 0.0), nrow=length(d0)) %*% 
            t(Zsvd[["u"]])
        stopifnot( isSymmetric(Ltmp) )
        stopifnot( identical(dim(Z), dim(Ltmp)) )

        # `lambda2Set` does not contain the diagonal entries in `Z-Ltmp`
        lambda2Set <- sort(abs( (Z - Ltmp)[upper.tri(Z)] ), decreasing=TRUE)[1:maxCard]
        
        for (lambda2 in lambda2Set) {
            drrobj <- solve_square_drr(Z, lambda1, lambda2, S.diag.penalize, tol, maxIter)
            
            L0 <- drrobj[["L"]]
            S0 <- drrobj[["S"]]
            rankL <- sum(abs(drrobj[["Lsv"]]) >= valuetol)
            cardS <- sum(abs(drrobj[["S"]]) >= valuetol)

            drrErr <- (Z - L0 - S0)^2
            loss <- sum(drrErr[upper.tri(drrErr, diag=TRUE)])
            ebic <- -2*log(loss) + cardS*(log(N) + 2*ebic.gamma*log(p*(p+1)/2))
            isConvergent <- drrobj[["isConvergent"]]
            drrobj$loss <- loss / Z.norm

            performance[k, ] <- c(lambda1, lambda2, rankL, cardS, loss, ebic, isConvergent)
            k <- k + 1
        }
    }

    optResult <- as.data.table(performance) # [cardS <= maxCard]
    nearestRank <- optResult[which.min(abs(rankL-rstar)), rankL]
    optResult <- optResult[rankL==nearestRank][which.min(loss)]
    optDRR <- solve_square_drr(Z, optResult[["lambda1"]], optResult[["lambda2"]], S.diag.penalize, tol, maxIter)
    
    L0 <- drrobj[["L"]]
    S0 <- drrobj[["S"]]
    drrErr <- (Z - L0 - S0)^2
    optDRR$loss <- sum(drrErr[upper.tri(drrErr, diag=TRUE)]) / Z.norm

    stopifnot( identical(
        c(optResult[["lambda1"]], optResult[["lambda2"]]), 
        c(optDRR[["lambda1"]], optDRR[["lambda2"]]) 
    ) )

    attr(optDRR, "rankStar") <- rstar
    class(optDRR) <- "drr_factor_model"
    
    return(optDRR)
}
