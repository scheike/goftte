#ifndef CUMRES_H
#define CUMRES_H

#include <cmath>
#include <omp.h>
#include "matrix.h" 


inline double ss2(const scythe::Matrix<double> &M) {
  return(std::sqrt(scythe::sum(M%M)));
}

template <typename T>
inline scythe::Matrix<T> cumsum2(const scythe::Matrix<T> &M, const scythe::Matrix<T> &ord) {
  scythe::Matrix<T> res = cumsum(M);
  unsigned n = M.rows();  
  if (M.rows()==1)
    n = M.cols();    
  for (unsigned i=(n-1); i>0; i--) {
    if (fabs(ord[i]-ord[i-1])<1e-12) {
      res(i-1,scythe::_) = res(i,scythe::_);
    }
  }
  return(res);
}

double KolmogorovSmirnov(const scythe::Matrix<double> &W);

double CramerVonMises(const scythe::Matrix<double> &W, const scythe::Matrix<double> &It, unsigned var);

double AndersonDarling (const scythe::Matrix<double> &W, const scythe::Matrix<double> &It, unsigned var);

scythe::Matrix<double> ProdMat(const scythe::Matrix<double> &A, const scythe::Matrix<double> &B);

scythe::Matrix<double> Kz(const scythe::Matrix<double> &x, const double b);

scythe::Matrix<double> Wi(const scythe::Matrix<double> &r, 
			  const scythe::Matrix<double> &x,
			  const scythe::Matrix<double> &betaiid,
			  const scythe::Matrix<double> &etaraw, const double b, 
			  scythe::Matrix<double> &Wobs);

scythe::Matrix<double> Cpred(const scythe::Matrix<double> &cum, 
				    const scythe::Matrix<double> &x);
            
scythe::Matrix<> SumMat(const scythe::Matrix<> &M);

scythe::Matrix<> cumsum(const scythe::Matrix<> &M);

#endif /* CUMRES_H */
