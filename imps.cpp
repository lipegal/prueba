#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boos/numeric/ublas/io.hpp>
#include <lapack.h>
#include <matrixtypes.h>
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace ula;

int m=3; //Internal degrees of freedom (spin)
int M=7; //Virtual dimension
int L;  //Lattice sites
ComplexMatrix TrnsfrMtrx(ComplexMatrix &A, int k){
  int alpha,beta,Ns,alpha_,beta_;
  SubComplexMatrix Ak_(A,boost::numeric::ublas::range(0,m*M),boost::numeric::ublas::range(k*M,(k+1)*M));
  ComplexMatrix Ak(m,M*M);
  ComplexMatrix E_(M*M,M*M);
  ComplexMatrix E(M*M,M*M);
  for(Ns=0;Ns<m;Ns++)
    for(alpha=0;alpha<M;alpha++)
      for(beta=0;beta<M;beta++){
Ak(Ns,alpha*M+beta)=Ak_(Ns*M+alpha,beta);
      }
  noalias(E_)=prod(herm(Ak),Ak);
  for(alpha=0;alpha<M;alpha++)
    for(beta=0;beta<M;beta++)
      for(alpha_=0;alpha_<M;alpha_++)
for(beta_=0;beta_<M;beta_++)
 E(alpha*M+alpha_,beta*M+beta_)=E_(alpha*M+beta,alpha_*M+beta_);
  return E;
  }

ComplexMatrix TrnsfrMtrx_(ComplexMatrix &A,ComplexMatrix &B, int k){
  int alpha,beta,Ns,alpha_,beta_;
  SubComplexMatrix Ak_(A,boost::numeric::ublas::range(0,(m+1)*M),boost::numeric::ublas::range(k*M,(k+1)*M));
  SubComplexMatrix Bk_(B,boost::numeric::ublas::range(0,(m+1)*M),boost::numeric::ublas::range(k*M,(k+1)*M));
  ComplexMatrix Ak(m+1,M*M);
  ComplexMatrix Bk(m+1,M*M);
  ComplexMatrix E_(M*M,M*M);
  ComplexMatrix E(M*M,M*M);
  for(Ns=0;Ns<m;Ns++)
    for(alpha=0;alpha<M;alpha++)
      for(beta=0;beta<M;beta++){
Ak(Ns,alpha*M+beta)=Ak_(Ns*M+alpha,beta);
Bk(Ns,alpha*M+beta)=Bk_(Ns*M+alpha,beta);
      }
  noalias(E_)=prod(herm(Ak),Bk);
  for(alpha=0;alpha<M;alpha++)
    for(beta=0;beta<M;beta++)
      for(alpha_=0;alpha_<M;alpha_++)
for(beta_=0;beta_<M;beta_++)
 E(alpha*M+alpha_,beta*M+beta_)=E_(alpha*M+beta,alpha_*M+beta_);
  return E;
}
int main(){


return 0;
}
