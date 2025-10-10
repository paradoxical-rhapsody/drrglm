#' @name nrr
#' @title Nuclear-Norm Regularized Matrix-Variate GLMs
#' @description 
#' Solve the following problem
#' \deqn{
#'  \min_{L} \Big\lbrace 
#'      \frac{1}{2N} \sum_{n=1}^N \ell(y_n, X_n) + \lambda \|L\|_* 
#'  \Big\rbrace.
#' }
#' 
#' @param x \eqn{p_1 \times p_2 \times N} numeric array with **mean zero**.
#' @param y Numeric vector of length \eqn{N}.
#' @param family See [glm] and [family].
#' @param lambda \eqn{\lambda}.
#' @param L0 Initialization to \eqn{L}.
#' @param tol Convergence tolerance.
#' @param maxIter Maximal step of iterations.
#' @param verbose Print iterations?
#' 
#' @return `list(L, Lsv, intercept, iter, isConvergent, lambda, deltaval)`.
#' 
#' @examples 
#' set.seed(2025)
#' r0 <- 13
#' c0 <- 17
#' N <- 500
#' x <- array(runif(r0*c0*N), c(r0, c0, N))
#' y <- rnorm(N)
#' 
#' family <- "gaussian"
#' lambda <- 0.15
#' 
#' system.time( result <- nrrglm(x, y, family, lambda, verbose=TRUE) )
#' 
#' @noRd 
nrrglm <- function(x, y, family, lambda, 
                L0=NULL, tol=1e-3, maxIter=500, verbose=FALSE){
    stopifnot( length(dim(x)) == 3 )
    stopifnot( length(y) == dim(x)[3] )
    stopifnot( lambda > 0 )
    stopifnot( length(lambda) == 1 )

    N <- length(y)
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

    # fml <- formula(y ~ 0)
    # if (intercept) fml <- formula(y ~ 1)
    loss <- function(C) {
        eta <- apply(x, 3, function(xn) sum(C * xn))
        # -1.0/N * logLik(glm(y ~ 0, family=family, offset=eta))
        1.0/(2*N) * deviance(glm(y ~ 0, family=family, offset=eta))
    }

    grad <- function(C) {
        eta <- apply(x, 3, function(xn) sum(C * xn))
        mu <- family$linkinv(eta)
        w <- (mu - y) * family$mu.eta(eta) / family$variance(mu)
        apply( sweep(x, 3, w, "*"), 1:2, mean )
    }

    showiter <- function(verbose)
        if (verbose) message(sprintf("iter %i: obj=%.4f", k, obj))

    if (is.null(L0))
        L0 <- matrix(0.0, NROW(x), NCOL(x))

    W <- L <- L0
    L.sv <- svd(L)
    svals <- L.sv[["d"]]
    obj <- obj0 <- loss(L) + lambda*sum(svals) 
    
    alpha0 <- 1.0
    alpha <- (1.0 + sqrt(1.0 + 4*alpha0^2)) / 2.0
    iter <- rep(NA_real_, maxIter)
    isConvergent <- FALSE
    delta.factor <- 1.0
    for (k in seq_len(maxIter)) {
        showiter(verbose)
        iter[k] <- obj

        ## --- SEARCH POINT ---
        W <- L + (alpha0-1.0)/alpha * (L - L0)
        deltaval <- delta.factor

        ## --- APPROXIMATE SOLUTION ---
        L0 <- L
        obj0 <- obj
        W.go <- grad(W) / deltaval
        if ( max(abs(W.go)) < min(1.0e-8, 0.1*tol) ) break

        L.sv <- svd(W - W.go)
        svals <- pmax(L.sv[["d"]] - lambda/deltaval, 0.0)
        idx <- union( which.max(svals), which(svals > 0.0) )
        L <- L.sv[["u"]][, idx] %*% diag(svals[idx], nrow=length(idx)) %*% t(L.sv[["v"]][, idx])
        obj <- loss(L) + lambda*sum(svals)
        stopifnot( svals >= 0 )
        stopifnot( length(idx) >= 1 )
        stopifnot( length(obj) == 1 )
        stopifnot( is.numeric(obj) )

        ## --- FORCE DESCENT ---
        if (obj > obj0) {
            L <- L0
            obj <- obj0
            delta.factor <- 1.5 * delta.factor
            next
        }

        ## --- CHECK CONVERGENCE ---
        if ( abs(obj-obj0) <= (abs(obj0) + 0.1)*tol ) {
            isConvergent <- TRUE
            break
        }

        ## --- UPDATE PARAs ---
        alpha0 <- alpha
        alpha <- (1.0 + sqrt(1.0 + 4*alpha0^2)) / 2.0
    }

    # coef.intercept <- 0.0
    # if (intercept) {
    #     eta <- apply(x, 3, function(xn) sum(L * xn))
    #     coef.intercept <- coef(glm(fml, family=family, offset=offset))
    # }
    # stopifnot( length(coef.intercept) == 1 )

    egg <- list(
        L = L, 
        Lsv = svals,
        # intercept = coef.intercept, 
        iter = c(na.omit(iter)),
        isConvergent = isConvergent,
        lambda = lambda,
        deltaval = deltaval
    )
    class(egg) <- "nrrglm"

    return(egg)
}
