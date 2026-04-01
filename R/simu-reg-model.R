#' @name simu-reg-model
#' @title Simulation: Matrix-Variate GLM
#' 
#' @param N Sample size.
#' @param p1 Row dimension.
#' @param p2 Column dimension.
#' @param C0 Coefficient matrix.
#' @param family See [glm] and [family].
#' @param rank0 Rank of `L`.
#' @param S.type Type of `S`.
#' \itemize{
#'  \item `zero`: Zero matrix.
#'  \item `rand`: Nonzero elements at random positions.
#'  \item `block`: Block matrix.
#'  \item `diag`: Diagonal matrix.
#' }
#' @param err.student.dof The degree of freedom of Student-\eqn{t} for error term in
#'  `family="gaussian"`. By default it is `NULL` for normal error.
#' @param seed Random seed.
#' 
#' @return
#' \itemize{
#'  \item `simu_reg_coefs`: `list(L, S)`.
#'  \item `simu_reg_data`: `list( train=list(x, y), test=list(x, y) )`.
#' }
#' 
#' @examples
#' p1 <- 10
#' p2 <- 9
#' rank0 <- 3
#' 
#' coefs <- simu_reg_coefs(p1, p2, rank0, "zero")
#' coefs <- simu_reg_coefs(p1, p2, rank0, "rand")
#' coefs <- simu_reg_coefs(p1, p2, rank0, "block")
#' coefs <- simu_reg_coefs(p1, p2, rank0, "diag")
#' 
#' N <- 100
#' S.type <- "rand"
#' family <- "binomial"
#' 
#' coefs <- simu_reg_coefs(p1, p2, rank0, S.type)
#' C0 <- coefs[["L"]] + coefs[["S"]]
#' 
#' set.seed(2025)
#' DT <- simu_reg_data(family, N, C0)
#' DT <- simu_reg_data("gaussian", N, C0, err.student.dof=3)
#' 
NULL



#' @rdname simu-reg-model
#' @order 1
#' @export
simu_reg_coefs <- function(p1, p2, rank0, S.type, seed=2023){
    if (!is.null(seed)) set.seed(seed)

    stopifnot( S.type %in% c("zero", "rand", "block", "diag") )
    stopifnot( rank0 >= 0 && rank0 <= min(p1, p2) )
    stopifnot( min(p1, p2) >= 4 )
    stopifnot( p1*p2 >= 10 )
    
    L <- matrix(rnorm(p1*p2), p1)
    Lsvd <- svd(L)
    u <- Lsvd[["u"]]
    v <- Lsvd[["v"]]
    sval <- c( runif(rank0, 8.0, 10.0), rep(0.0, min(p1, p2)-rank0) )
    L <- u %*% diag(sval, nrow=min(p1, p2)) %*% t(v)

    S <- matrix(0.0, p1, p2)
    idx <- switch(S.type,
        "zero" = NULL,
        "rand" = sample(p1*p2, 10),
        "block" = which( row(S) <= 4 & col(S) <= 4 ), 
        "diag" = which(diag(TRUE, p1, p2))
    )
    S[idx] <- runif(length(idx), 8.0, 10.0)

    return( list(L=L, S=S) )
}



#' @rdname simu-reg-model
#' @order 2
#' @export
simu_reg_data <- function(family, N, C0, err.student.dof=NULL){
    stopifnot( family %in% c('gaussian', 'binomial') )
    rDim <- NROW(C0)
    cDim <- NCOL(C0)

    set_x_and_y <- function(family, N){
        x <- array(rnorm(rDim*cDim*N), c(rDim, cDim, N))
        eta <- apply(x, 3, function(xn) sum(xn*C0))

        if (family == "gaussian") {
            if (is.null(err.student.dof)) {
                y <- eta + rnorm(N, sd=sd(eta)/8)
            } else {
                y <- eta + rt(N, err.student.dof)
            }
        }
        if (family == "binomial") {
            y <- rbinom(N, 1, binomial()$linkinv(eta))
        }

        stopifnot( length(y) == dim(x)[3] )
        return( list(x=x, y=y) )
    }
    
    list(
        train = set_x_and_y(family, N),
        test  = set_x_and_y(family, 1000)
    )
}
