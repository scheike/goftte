#ifndef EXTRA_H
#define EXTRA_H

#include <cmath>
#include <functional>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>
//#include <scythestat/rng/mersenne.h> 
//#include <scythestat/distributions.h> 
//#include <scythestat/ide.h> 
//#include <scythestat/la.h> 
#include "matrix.h"
//#include <scythestat/rng.h> 
//#include <scythestat/smath.h> 
//#include <scythestat/stat.h> 
//#include <scythestat/optimize.h> 


using namespace scythe;
using namespace std;

/* inline scythe::Matrix<double> loadMatrix(char filename[]) { */
/*   scythe::Matrix<double> matrix; */
/*   std::ifstream in(filename); */
/*   in >> matrix;   */
/*   //  return old; */
/* } */


/* WE NEED A FUNCTION WHICH RETURNS A UNIQUE ORDERING (order RETURNS RANKS! e.g.
   4 1 1 1 5 instead of 4 1 2 3 5 */
/* 
template <typename T, matrix_order O, matrix_style S>
Matrix<unsigned int, O, Concrete> myorder (const Matrix<T,O,S>& M)
{
   Matrix<unsigned int, O, Concrete> res = order<O,Concrete>(M);
  //  Matrix<unsigned int, O, Concrete> res = order_alg<O,Concrete>(M);
  unsigned prev = res[0];
  for (unsigned i=1; i<res.size(); i++) {
    unsigned cur = res[i];
    //    if (prev=cur)
    prev = cur;
  }
  return(res);
}
*/

inline unsigned sum(const Matrix<bool> &m) {
  unsigned count=0;
  for (unsigned i=0; i<m.size(); i++)
    if(m[i]) count++;
  return(count);
}
inline Matrix<unsigned> boolidx(const Matrix<bool> &idx) {
  unsigned n = sum(idx);
  Matrix<unsigned> Ridx(n,1);
  unsigned pos = 0;
  for (unsigned i=0; i<n; i++)  {
    if (idx[i]) { Ridx[i]=pos; pos++; }
  }
  return(Ridx);
}

template <typename T>
inline Matrix<T> chrows(const Matrix<T> &m, const Matrix<bool> &idx) {
  assert (m.rows()==idx.size()); 
  Matrix<unsigned> Ridx = boolidx(idx);
  Matrix<T> res(Ridx.size(), m.cols());
  for (unsigned i=0; i<Ridx.size(); i++) {
    res(i,_) = m(Ridx[i],_);
  }  
  return(res);
}
template <typename T>
inline Matrix<T> chrows(const Matrix<T> &m, const Matrix<unsigned> &idx) {
  Matrix<T> res(idx.size(), m.cols(), false);
  for (unsigned i=0; i<idx.size(); i++) {
    res(i,_) = m(idx[i],_);
  }
  return(res);
}
template <typename T>
inline Matrix<T> chcols(const Matrix<T> &m, const Matrix<bool> &idx) {
  assert (m.cols()==idx.size()); 
  std::vector<unsigned> Ridx;
  for (unsigned i=0; i<m.cols(); i++)  {
    if (idx[i]) { Ridx.push_back(i); }
  }
  Matrix<T> res(Ridx.size(), m.cols());
  for (unsigned i=0; i<Ridx.size(); i++) {
    res(_,i) = m(_,Ridx[i]);
  }  
  return(res);
}
template <typename T>
inline Matrix<T> chcols(const Matrix<T> &m, const Matrix<int> &idx) {
  Matrix<T> res(m.rows(), idx.size(), false);
  for (unsigned i=0; i<idx.size(); i++) {
    res(_,i) = m(_,idx[i]);
  }
  return(res);
}

template <typename T>
inline Matrix<T> addCol(const Matrix<T> &m, const Matrix<T> &c) { // Add vector c to each column of m
  Matrix<T> res = m;
  for (unsigned i=0; i<m.cols(); i++)
    res(_,i) += c;
  return(res);
}
template <typename T>
Matrix<T> addRow(const Matrix<T> &m, const Matrix<T> &c) { // Add vector c to each column of m
  Matrix<T> res = m;
  for (unsigned i=0; i<m.rows(); i++)
    res(i,_) += c;
  return(res);
}
template <typename T>
Matrix<T> multCol(const Matrix<T> &m, const Matrix<T> &c) { // Multiply vector c with each column of m
  Matrix<T> res = m;
  for (unsigned i=0; i<m.cols(); i++)
    res(_,i) %= c;
  return(res);
}
template <typename T>
inline Matrix<T> multRow(const Matrix<T> &m, const Matrix<T> &c) { // Multiply vector c with each column of m
  Matrix<T> res = m;
  for (unsigned i=0; i<m.rows(); i++)
    res(i,_) %= c;
  return(res);
}


template <typename T>
inline Matrix<T> reverse(const Matrix<T> &m) {
  Matrix<T> res = m;
  reverse(res.begin(), res.end());
  return(res);
}

template <typename T>
inline Matrix<T> sapply(const Matrix<T> &m, T (*fun)(T)){
  Matrix<T> res = m;
  for (unsigned i=0; i<m.rows(); i++) {
    for (unsigned j=0; j<m.cols(); j++) {
      res(i,j) = (*fun)(m(i,j));      
    }
  } 
  return(res);  
}


template <typename T>
inline Matrix<T> apply2(const Matrix<T> &m, unsigned doRow, Matrix<T> (*fun)(const Matrix<T>&)){
  if (doRow==1) {
    Matrix<T> r0 = (*fun)(m(0,_));
    return(r0);  
  } else {
    Matrix<T> r0 = (*fun)(m(_,0));
    return(r0);  
  }
}

//    cerr << apply(X,1, scythe::sum) << endl;
//    cerr << apply(X,2, scythe::sum) << endl;
template <typename T>
inline  Matrix<T> apply(const Matrix<T> &m, unsigned doRow,T (*fun)(const Matrix<T>&)){
  if(doRow==1){
    Matrix<T> res(m.rows(),1);
    for (unsigned i=0; i<m.rows(); i++) {
      Matrix<T> r = m(i,_);
      res[i] = (*fun)(r);      
      //      res(i,_) = fun(r);
    } 
    return(res);  
  } else {
    Matrix<T> res(m.cols(),1);
    for (unsigned i=0; i<m.cols(); i++) {
      Matrix<T> r = m(_,i);
      //      res(i,_) = fun(r);
      res[i] = (*fun)(r);
    }
    return(res);   
  }
}

template <typename T, typename FUNCTOR>
inline Matrix<T> apply(const Matrix<T> &m, unsigned doRow, FUNCTOR fun){
  //Matrix<T> apply(const Matrix<T> &m, unsigned doRow, T (*fun)(const Matrix<T>&)){
  if(doRow==1){
    Matrix<T> res(m.rows(),1);
    for (unsigned i=0; i<m.rows(); i++) {
      Matrix<T> r = m(i,_);
      //      res[i] = (*fun)(r);      
      res[i] = fun(r);
    } 
    return(res);  
  } else {
    Matrix<T> res(m.cols(),1);
    for (unsigned i=0; i<m.cols(); i++) {
      Matrix<T> r = m(_,i);
      res[i] = fun(r);
    }
    return(res);   
  }
}

 
template <typename T>
inline std::string Rout2(scythe::Matrix<T> const& M) {
  std::ostringstream out;
  
  for(unsigned r = 0; r < M.rows(); ++r) {
    out << std::endl << "[" << r+1 << "] ";
    //    out << "";      
    for(unsigned int s = 0; s < M.cols(); ++s) {
      if(s > 0) out << ", ";
      //      out << (*this)(r,s);
      out << M(r,s);
    }
  }
  out << std::endl;
  return std::string(out.str());
}

template <typename T>
inline std::string Rout(scythe::Matrix<T> const& M) {
  std::ostringstream out;
  out << "rbind("; //<< std::endl;
  for(unsigned r = 0; r < M.rows(); ++r) {
    if(r > 0) out << "), " << std::endl;
    out << "c(";      
    for(unsigned int s = 0; s < M.cols(); ++s) {
      if(s > 0) out << ", ";
      //      out << (*this)(r,s);
      out << M(r,s);
    }
  }
  out << "))" << std::endl;
  return std::string(out.str());
}

template <typename T>
inline std::string Rout(std::vector<T> const& x) {
  std::ostringstream out;
  out << "c(";
  for(unsigned r = 0; r < x.size(); ++r) {
    if(r > 0) out << ", ";
    //      out << (*this)(r,s);
    out << x[r];
  }
  out << ")" << std::endl;
  return std::string(out.str());
}


#endif /* EXTRA_H */
