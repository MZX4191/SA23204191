//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
//' @title A Gibbs sampler for multivariate normal distribution using Rcpp
//' @description A Gibbs sampler for multivariate normal distribution using Rcpp
//' @param m the number of samples
//' @param x0 the initial value
//' @param mu the mean vector of the distribution
//' @param sigma the covariance matrix of the distribution
//' @return a MCMC chain \code{n}
//' @examples
//' \dontrun{
//' a<-rmnGibbsc(2, c(1,2,3),c(0,0,0),matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),3,3))
//' a
//' }
//' @export
//[[Rcpp::export]]

arma::mat rmnGibbsc(int m, arma::vec x0, arma::vec mu,arma::mat sigma){
  int n=mu.size();
  arma::mat x(n,m) ;
  x.fill(0);  
    x.col(0)=x0;
    for(int i=1;i<m;i++){
      x.col(i)=x.col(i-1);
      for(int j=0 ;j<=(n-1);j++){
        arma::mat A(n,n);
        A.fill(0);
        for(int k=0;k<=(n-1);k++){
          for(int t=0;t<=(n-1);t++){
            if(k==t&&k<=(j-1)&&t<=(j-1)) A(k,t)=1;
            if(t>=(j+1)&&t==(k+1)) A(k,t)=1;
            if(k==(n-1)&&t==j) A(k,t)=1;
          }
        }
        arma::vec Amu=A*mu,Ax=A*x.col(i);
        arma::mat Sigma=A*sigma*A.t();
        arma::mat Sigma11=Sigma(arma::span(0,n-2),arma::span(0,n-2));
        arma::mat Sigma21=Sigma(n-1,arma::span(0,n-2));
        arma::mat Sigma12=Sigma(arma::span(0,n-2),n-1);
        arma::mat a=Amu(n-1)+Sigma21*arma::inv(Sigma11)*(Ax(arma::span(0,n-2))-Amu(arma::span(0,n-2)));
        arma::mat b=Sigma(n-1,n-1)-Sigma21*arma::inv(Sigma11)*Sigma12;
        Rcpp::NumericVector c=Rcpp::rnorm(1,a(0,0),b(0,0));
        x(j,i)=c(0);
      }
    }
    return(x);
}