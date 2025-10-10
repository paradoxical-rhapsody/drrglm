#' @name simu-factor-model
#' @title Simulation: Factor Model with Correlated Noise
#' 
#' @description 
#' \eqn{y = B u + w} with \eqn{u \sim N(0, I_{r})}, 
#' \eqn{w \sim N(0, S)} and \eqn{B \in \mathbb{R}^{p \times r}}.
#' 
#' \itemize{
#'  \item `simu_factor_model_paras`: simulate `(B, S)`.
#'  \item `simu_factor_model_data`: simulate `y`.
#' }
#' 
#' @param N Sample size.
#' @param p Dimension of data.
#' @param r Number of factors.
#' @param noise Mode of \eqn{S}:
#' \itemize{
#'  \item `diag`: Diagonal matrix.
#'  \item `rand`: Nonzero elements at random positions.
#'  \item `tridiag`: Tri-diagonal matrix.
#'  \item `block`: Block matrix.
#' }
#' @param B Matrix \eqn{p \times r}.
#' @param S Symmetric matrix \eqn{p \times p}.
#' 
#' @return 
#' \itemize{
#'  \item `simu_factor_model_paras`: `list(B, S)`.
#'  \item `simu_factor_model_data`: Numeric matrix `y` of \eqn{N \times p}.
#' }
#' 
#' 
#' @examples 
#' set.seed(2025)
#' N <- 50
#' p <- 15
#' r <- 7
#' 
#' paras <- simu_factor_model_paras(p, r, 'diag')
#' paras <- simu_factor_model_paras(p, r, 'rand')
#' paras <- simu_factor_model_paras(p, r, 'tridiag')
#' paras <- simu_factor_model_paras(p, r, 'block')
#' 
#' DT <- simu_factor_model_data(N, paras[["B"]], paras[["S"]])
#' 
NULL


#' @rdname simu-factor-model
#' @order 1
#' @export 
simu_factor_model_paras <- function(p, r, noise){
    stopifnot( p > r )
    stopifnot( noise %in% c('diag', 'rand', 'tridiag', 'block') )

    set.seed(2023)

    tmp <- svd( matrix(rnorm(p*r), p) )
    U <- tmp[["u"]][, 1:r, drop=FALSE]
    V <- tmp[["v"]]
    stopifnot( dim(U) == c(p, r) )
    stopifnot( dim(V) == c(r, r) )
    
    B <- U %*% diag(runif(r, 8, 10), nrow=r) %*% t(V)
    stopifnot( dim(B) == c(p, r) )

    idx <- switch(noise,
        diag = diag(TRUE, p), 
        rand = { 
            tmpidx <- diag(TRUE, p)
            tmpidx[sample(which(lower.tri(tmpidx)), 10)] <- TRUE
            tmpidx | t(tmpidx)
        }, 
        tridiag = { abs(row(diag(p)) - col(diag(p))) <= 1 }, 
        block = { pmax(row(diag(p)), col(diag(p))) <= 4 }
    )
    stopifnot( isSymmetric(idx) )
    switch(noise,
        diag = stopifnot( sum(idx) == p ),
        rand = stopifnot( sum(idx) == (p+20) ),
        tridiag = stopifnot( sum(idx) == (p + 2*(p-1)) ),
        block = stopifnot( sum(idx) == 16 )
    )
    diag(idx) <- FALSE

    E <- idx * matrix(runif(p*p, 0.5, 1.0), p, p)
    E[upper.tri(E, diag=TRUE)] <- 0.0
    E <- E + t(E)
    stopifnot( diag(E) == 0.0 )
    stopifnot( isSymmetric(E) )
    stopifnot( setequal( which(E != 0), which(idx) ) )

    c0 <- 0.9 * max(abs( B %*% t(B) ))
    if (noise == "diag") {
        stopifnot( all(E == 0) )
        S <- diag(c0, nrow=p)
    } else {
        stopifnot( any(E !=0) )
        S <- diag(1.1*abs(min(eigen(E)$values)), nrow=p) + E
        S <- (c0 / max(abs(S))) * S
    }
    stopifnot( isSymmetric(S) )

    return( list(B=B, S=S) )
}


#' @rdname simu-factor-model
#' @order 2
#' @export
simu_factor_model_data <- function(N, B, S){
    stopifnot( isSymmetric(S) )
    stopifnot( NROW(B) == NROW(S) )

    p <- NROW(B)
    r <- NCOL(B)
    y <- (matrix(rnorm(N*r), N) %*% t(B)) + (matrix(rnorm(N*p), N) %*% chol(S))
    stopifnot( dim(y) == c(N, p) )

    return(y)
}
