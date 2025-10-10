#' @title Initialize Parameters
#' @description Initialize \eqn{\widehat{L}} via nuclear-norm 
#' regularized regression.
#' 
#' @param x Numeric array.
#' @param y Numeric vector.
#' @param family See [glm] and [family].
#' @param maxRank The maximum of rank to be detected.
#' @param tol Convergence tolerance.
#' @param verbose Print iterations?
#' 
#' @return `list(
#'  family, 
#'  grad, 
#'  lipschitz, 
#'  sigma0, 
#'  rankStar, 
#'  rankStarHat, 
#'  lambda.star, 
#'  L0.star
#' )`
#' 
#' 
#' @export
ini_paras <- function(x, y, family, 
                    maxRank=min(floor(NROW(x)/2), floor(NCOL(x)/2)), 
                    tol=1e-3, verbose=FALSE){
    stopifnot( length(dim(x)) == 3 )
    stopifnot( length(y) == dim(x)[3] )

    valuetol <- min(0.1*tol, 1e-4) # 1e-6
    maxIter <- 500

    x.norm2 <- apply(x, 3, norm, type="F")^2
    stopifnot( length(x.norm2) == length(y) )

    if (is.character(family))
		family <- get(family, mode = "function", envir = parent.frame())
	if (is.function(family))
		family <- family()
	if (is.null(family$family)) {
		print(family)
		stop("'family' not recognized")
	}

    grad <- function(L) {
        eta <- apply(x, 3, function(xn) sum(L * xn))
        mu <- family$linkinv(eta)
        w <- (mu - y) * family$mu.eta(eta) / family$variance(mu)
        stopifnot( length(w) == length(eta) )

        apply( sweep(x, 3, w, "*"), 1:2, mean )
    }
    lipschitz <- function(L) { # Lipschitz const of gradient
        eta <- apply(x, 3, function(xn) sum(L * xn))
        mu <- family$linkinv(eta)
        stopifnot( length(eta) == length(mu) )
        stopifnot( length(eta) == length(x.norm2) )

        mean( (family$mu.eta(eta))^2 / family$variance(mu) * x.norm2 )
    }

    ## --- INITIALIZATION ---
    # L0 <- apply(x, 1:2, function(xij) mean(xij*y) / mean(xij^2))
    L0 <- matrix(0.0, NROW(x), NCOL(x))
    tmp <- nrrglm(x, y, family, 1e-4, L0, 0.1*tol, maxIter, FALSE)
    L0 <- tmp[["L"]]
    stopifnot( dim(L0) == c(NROW(x), NCOL(x)) )
    
    ## extend `sigma0` of length-`min` to length-`max`, used
    ##  in `dof` of OLS (Zhouetal2014)
    sigma0 <- tmp[["Lsv"]]
    if (length(sigma0) < max(NROW(x), NCOL(x))) {
        sigma0[max(NROW(x), NCOL(x))] <- NA_real_
        sigma0[is.na(sigma0)] <- 0.0
    }
    stopifnot( length(sigma0) == max(NROW(x), NCOL(x)) )

    c0 <- quantile(tmp[["Lsv"]], 0.10)
    rankStar <- which.max( get_gr_rank(tmp[["Lsv"]], c0) )
    rankStar <- min(rankStar, maxRank)
    stopifnot( is.numeric(rankStar) )
    stopifnot( length(rankStar) == 1 )

    lipschitzVal <- lipschitz(L0)
    lambda <- lipschitzVal * svd(L0 - grad(L0)/lipschitzVal)$d[rankStar]
    stopifnot( lambda > 0 )
    stopifnot( length(lambda) == 1 )

    ## --- SEARCH `lambda1Star` For `rankStar` ---
    lambdaLower <- (-Inf)
    lambdaUpper <- (+Inf)
    k <- 0
    while (TRUE) {
        k <- k + 1
        tmp <- nrrglm(x, y, family, lambda, L0, tol, maxIter, verbose=FALSE)
        L0 <- get("L", tmp)
        stopifnot( dim(L0) == c(NROW(x), NCOL(x)) )

        rankL <- sum(abs(tmp[["Lsv"]]) >= valuetol)
        
        if (verbose) 
            message(sprintf("%i-th: [%.6f, %.6f], rankL=%i(%i)", k, lambdaLower, lambdaUpper, rankL, rankStar))

        if (rankL == rankStar) break
        
        if (rankL < rankStar) lambdaUpper <- min(lambda, lambdaUpper)
        if (rankL > rankStar) lambdaLower <- max(lambda, lambdaLower)

        if ( abs(lambdaUpper - lambdaLower) < 1e-4 ) break

        stopifnot( lambdaLower < lambdaUpper )
        stopifnot( !is.na(lambdaLower) && !is.na(lambdaUpper) )
        stopifnot( is.finite(lambdaLower) || is.finite(lambdaUpper) )

        flag <- FALSE
        if (!is.finite(lambdaLower) && is.finite(lambdaUpper)) { # (-Inf, b)
            lambda <- 0.1 * lambdaUpper
            flag <- TRUE
        }
        if (is.finite(lambdaLower) && !is.finite(lambdaUpper)) { # (a, +Inf)
            lambda <- 10.0 * lambdaLower
            flag <- TRUE
        }
        if (is.finite(lambdaLower) && is.finite(lambdaUpper)) {  # (a, b)
            lambda <- (lambdaLower + lambdaUpper) / 2.0
            flag <- TRUE
        }
        stopifnot(flag)

        if (lambda < 0.01 * tol) break
        if (abs(lambdaUpper - lambdaLower) <= 0.1*tol) break
    }

    egg <- list(
        family = family,
        grad = grad,
        lipschitz = lipschitz,
        sigma0 = sigma0,
        rankStar = rankStar,
        rankStarHat = rankL,
        lambda.star = lambda,
        L0.star = get("L", tmp)
    )
    class(egg) <- "drr.ini"

    return(egg)
}
