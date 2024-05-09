#include <RcppEigen.h>
#ifdef _OPENMP
# include <omp.h>
#endif
using namespace Rcpp;
using namespace std;

//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::export]]
NumericMatrix gen_cost(NumericMatrix AR, NumericMatrix BR, int threads) {
  Eigen::Map<Eigen::MatrixXd> A(as<Eigen::Map<Eigen::MatrixXd> >(BR));
  Eigen::Map<Eigen::MatrixXd> B(as<Eigen::Map<Eigen::MatrixXd> >(AR));
  Eigen::setNbThreads(threads);
  int Al=A.rows();
  int Bl=B.rows();
  
  Eigen::MatrixXd x=(A*A.transpose()).diagonal();
  Eigen::MatrixXd z=(B*B.transpose()).diagonal();
  Eigen::MatrixXd onesA = Eigen::MatrixXd::Constant(1,Al, 1.0);
  Eigen::MatrixXd onesB = Eigen::MatrixXd::Constant(Bl,1, 1.0);
  Eigen::MatrixXd Cmat=(z*onesA+onesB*x.transpose()-(2.0*B*A.transpose()));
    // this can lead to (small) numerically negative values (and then NaNs when taking the root)
    
  // quick-and-(maybe dirty)-fix
  for(auto x : Cmat.reshaped()) {
    x = max(x, 0.0);
  }

  return  wrap(Cmat);
}


// naive alternative which is easier and about two times as fast
// (unless maybe if we use huge numbers of points and several threads)

//[[Rcpp::export]]
NumericMatrix gen_cost0(NumericMatrix xx, NumericMatrix yy) {
  int m = xx.nrow();
  int n = yy.nrow();
  
  NumericMatrix Cmat(m,n);  // initializes with zeros
  for (int j=0; j<n; j++) {
  for (int i=0; i<m; i++) {
    double d1 = xx(i,0) - yy(j,0);
    double d2 = xx(i,1) - yy(j,1);
    Cmat(i,j) = d1 * d1 + d2 * d2;
  }
  }
  
  return Cmat;
}


//[[Rcpp::export]]
NumericMatrix gen_cost0d(NumericMatrix xx, NumericMatrix yy) {
  int m = xx.nrow();
  int n = yy.nrow();
  int d = xx.ncol();
  if (d != yy.ncol()) {
    stop("number of coordinates must agree");
  }
  
  NumericVector temp(d);
  NumericMatrix Cmat(m,n);  // initializes with zeros
  for (int j=0; j<n; j++) {
    for (int i=0; i<m; i++) {
      for (int k=0; k<d; k++) {
        temp[k] = xx(i,k) - yy(j,k);
        Cmat(i,j) += temp[k] * temp[k];
      }
    }
  }
  
  return Cmat;
}
