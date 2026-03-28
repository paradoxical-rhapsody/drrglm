#' @name drr
#' @title Doubly Regularized Matrix-Variate GLMs
#' @description
#' Solve the following problem
#' \deqn{
#'   \min_{L, S} \Big\lbrace \frac{1}{N}
#'      \sum_{n=1}^N \ell(y_n, X_n) + \lambda_1 \|L\|_* + \lambda_2 \|S\|_1
#'   \Big\rbrace,
#' }
#' where \eqn{\ell(y_n, X_n)} refers to the negative-log-likelihood on \eqn{(y_n, X_n)} under pre-specified GLM.
#' 
#' In linear model, the problem has the form
#' \deqn{
#'      \min_{L, S} \Big\lbrace \frac{1}{2N} \sum_{n=1}^N
#'          \big( y_n - \text{tr}(X_n'C) \big)^2 +
#'          \lambda_1 \|L\|_* + \lambda_2 \|S\|_1 \Big\rbrace,
#'      ~~\text{s.t.}~~ C=L+S.
#' }
#' 
#' @param x \eqn{p_1 \times p_2 \times N} numeric array with **mean zero**.
#' @param y Numeric vector of length \eqn{N}.
#' @param family See [glm] and [family].
#' @param lambda1 \eqn{\lambda_1}.
#' @param lambda2 \eqn{\lambda_2}.
#' @param C0 Initialization to \eqn{C}.
#' @param tol Convergence tolerance.
#' @param maxIter Maximal step of iterations.
#' @param verbose Print iterations?
#' 
#' @return `list(L, S, Lsv, iter, lambda1, lambda2, isConvergent)`.
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
#' lambda1 <- 0.15
#' lambda2 <- 0.04
#' 
#' system.time( egg <- drrglm(x, y, family, lambda1, lambda2, verbose=TRUE) )
#' 
#' @export
drrglm <- function(x, y, family, lambda1, lambda2,
                C0=NULL, tol=1e-3, maxIter=300, verbose=FALSE){
    stopifnot( length(dim(x)) == 3 )
    stopifnot( length(y) == dim(x)[3] )
    stopifnot( lambda1 > 0 && lambda2 > 0 )
    stopifnot( length(lambda1) == 1 )
    stopifnot( length(lambda2) == 1 )

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

    if (is.null(C0))
        C0 <- matrix(0.0, NROW(x), NCOL(x))

    W <- C <- C0
    obj <- obj0 <- loss(C) + lambda2*sum(abs(C))

    alpha <- alpha0 <- 1.0
    iter <- rep(NA_real_, maxIter)
    isConvergent <- FALSE
    delta.factor <- 1.0
    for (k in seq_len(maxIter)) {
        showiter(verbose)
        iter[k] <- obj

        W <- C + (alpha0-1.0)/alpha * (C - C0)
        deltaval <- delta.factor

        C0 <- C
        obj0 <- obj
        W.go <- grad(W) / deltaval
        if ( max(abs(W.go)) < min(1.0e-8, 0.1*tol) ) break

        tmp <- cxx_solve_prox(W - W.go, lambda1/deltaval, lambda2/deltaval,
                            tol=0.5*tol, maxIter=1.5*maxIter)
        C <- tmp[["L"]] + tmp[["S"]]
        obj <- loss(C) + lambda1*sum(tmp[["Lsv"]]) + lambda2*sum(abs(tmp[["S"]]))
        stopifnot( length(obj) == 1 )
        stopifnot( is.numeric(obj) )

        if (obj > obj0) {
            C <- C0
            obj <- obj0
            delta.factor <- 1.5 * delta.factor
            next
        }

        if ( abs(obj-obj0) <= (abs(obj0) + 0.1)*tol ) {
            isConvergent <- TRUE
            break
        }

        alpha0 <- alpha
        alpha <- (1.0 + sqrt(1.0 + 4*alpha0^2)) / 2.0
    }

    egg <- list(
        L = tmp[["L"]],
        S = tmp[["S"]],
        Lsv = tmp[["Lsv"]],
        iter = c(na.omit(iter)),
        lambda1 = lambda1,
        lambda2 = lambda2,
        isConvergent = isConvergent
    )
    class(egg) <- "drrglm"

    return(egg)
}
