// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp ;


double cxx_huberloss(double x, double t) {
    if (fabs(x) <= t) {
        return x*x / 2.0 ;
    } else {
        return t * (fabs(x) - t/2.0) ;
    }
}



double cxx_hubergrad(double x, double t) {
    if (fabs(x) <= t) {
        return x ;
    } else {
        return copysign(t, x) ;
    }
}



double cxx_loss(arma::mat X, double t) {
    int rDim = X.n_rows ;
    int cDim = X.n_cols ;
    double egg = 0.0 ;
    for (int iR=0; iR < rDim; iR++) {
        for (int iC=0; iC < cDim; iC++) {
            egg += cxx_huberloss(X(iR, iC), t) ;
        }
    }
    return egg ;
}




arma::mat cxx_grad(arma::mat X, double t) {
    int rDim = X.n_rows ;
    int cDim = X.n_cols ;
    arma::mat egg(rDim, cDim) ;
    for (int iR=0; iR < rDim; iR++) {
        for (int iC=0; iC < cDim; iC++) {
            egg(iR, iC) = cxx_hubergrad(X(iR, iC), t) ;
        }
    }
    return egg ;
}




List cxx_soft_svd(arma::mat X, double t) {
    arma::mat u, v ;
    arma::vec svals ;
    
    svd(u, svals, v, X) ;
    svals = arma::clamp(svals - t, 0.0, arma::datum::inf) ;
    
    arma::mat L(size(X), arma::fill::zeros) ;
    // arma::uvec idx = find(svals > 0.0) ;
    for (unsigned int i=0; i < svals.n_elem; i++) {
        if (svals(i) > 0.0) {
            L += svals(i) * (u.col(i) * v.col(i).t()) ;
        }
    }
    
    return List::create(
        ::Named("L") = L, 
        ::Named("svals") = svals
    ) ;
}



// [[Rcpp::export]]
List cxx_solve_huber_nuclear(
                arma::mat Z, 
                double t, 
                double lambda, 
                double tol,
                int maxIter,
                bool verbose=false){
    arma::mat L=Z, L0=Z, W=Z ;
    arma::vec svals ;
    arma::mat Lu, Lv ;
    svd(Lu, svals, Lv, L) ;
    double obj = cxx_loss(L-Z, t) + lambda*sum(svals) ;

    double obj0 ;
    // arma::uvec idx ;
    // idx = find(svals > 0.0) ;

    double alpha0 = 1.0 ;
    double alpha = (1.0 + sqrt(1.0 + 4.0*alpha0*alpha0)) / 2.0 ;
    arma::vec iter(maxIter) ; 
    iter.fill(arma::datum::nan) ;
    bool isConvergent = false ;
    for (int k=0; k < maxIter; k++) {
        if (verbose) Rcout << "iter " << k+1 << " : obj=" << obj << "\n" ;
        iter(k) = obj ;

        // search point
        W = L + (alpha0 - 1.0)/alpha * (L - L0) ;

        // approximate solution
        L0 = L ;
        obj0 = obj ;

        svd(Lu, svals, Lv, W - cxx_grad(W-Z, t)) ;
        svals = arma::clamp(svals - lambda, 0.0, arma::datum::inf) ;
        L.fill(0.0) ;
        for (unsigned int i=0; i < svals.n_elem; i++) {
            if (svals(i) > 0.0) {
                L += svals(i) * (Lu.col(i) * Lv.col(i).t()) ;
            }
        }
        obj = cxx_loss(L-Z, t) + lambda*sum(svals) ;

        if (obj > obj0) {
            L = L0 ;
            obj = obj0 ;
            continue ;
        }

        if ( fabs(obj-obj0) <= (fabs(obj0)+0.1)*tol) {
            isConvergent = true ;
            break ;
        }

        // update paras
        alpha0 = alpha ;
        alpha = (1.0 + sqrt(1.0 + 4.0*alpha0*alpha0)) / 2.0 ;
    }
    
    arma::vec iter_good = iter.elem(find_finite(iter)) ;
    return List::create(
        ::Named("L") = L,
        ::Named("svals") = svals,
        ::Named("iter") = iter_good,
        ::Named("isConvergent") = isConvergent
    ) ;
}



arma::mat cxx_soft_threshold(arma::mat Z, double lambda){
    // legnth(lambda) == 1
    // lambda > 0
    arma::mat egg = sign(Z) % arma::max(arma::abs(Z) - lambda, arma::zeros<arma::mat>(size(Z))) ;
    return egg ;
}



// [[Rcpp::export]]
List cxx_solve_prox(arma::mat Z, double lambda1, double lambda2, 
    double tol, int maxIter){

    List Ltmp = cxx_solve_huber_nuclear(Z, lambda2, lambda1, tol, maxIter) ; 
    arma::mat S = cxx_soft_threshold(Z - as<arma::mat>(Ltmp["L"]), lambda2);

    return List::create(
        ::Named("L") = Ltmp["L"], 
        ::Named("S") = S,
        ::Named("Lsv") = Ltmp["svals"], 
        ::Named("iter") = Ltmp["iter"], 
        ::Named("isConvergent") = Ltmp["isConvergent"], 
        ::Named("lambda1") = lambda1, 
        ::Named("lambda2") = lambda2 
    ) ;
}
