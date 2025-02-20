/*

TeamMPS -Sebastian and Andres- Lib

This library contains many procedures needed for the time
evolution program for the ladder as well as for the chain.
It has to be included after the definition of the global program
parameters.

*/

#include <lapack.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <matrixtypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

using namespace ula;
using Range= boost::numeric::ublas::range; //G.

Complex I(0,1);  

// Kronecker delta
int delta(int i, int j){
  if(i==j)return 1;
  else return 0;
}

// For real numbers
int deltaR(double i, double j){
   if(fabs(i-j)<1e-7)return 1;
   else return 0;
}

//Factorial
int factorial(int num)
{
  if (num==0)
    return 1;
  if (num==1)
    return 1;
  return factorial(num-1)*num; // recursive call
}  

//Gaussian Function
Real Gauss(int x, int x0, Real sigma){
	return exp(-(x-x0)*(x-x0)/(2.*sigma*sigma))/(sigma*sqrt(2.*PI));
}

//-------------------------Operators------------------------//

//Identity operator
ComplexMatrix ident_op(){
  ComplexMatrix Id(m,m);
  for(int n=0;n<m;n++)
    Id(n,n)=1.;
  return Id;  
}


//Creation operator for different magnetic projections (sigma)
ComplexMatrix Creation(int sigma){
//Sigma {-1,0,1} -> {0,1,2}
	ComplexMatrix O(m,m);
	for(int i=0;i<m/4;i++){
		for(int j=0;j<m/4;j++){
			
			for(int alpha=0;alpha<m/4;alpha++){
				for(int beta=0;beta<m/4;beta++){
					int x = i*(m/4)+alpha;
					int y = j*(m/4)+beta;
					if(i==j){
						if(beta==0 && alpha>0 ){
							O(x,y)=delta(alpha-1,sigma);
							}
						}
					if(j==0 && i>0){
						if(alpha==beta){
							O(x,y)=delta(i-1,sigma);
							}
						}
					}
				}
			}	
		}
	return O;			
}

//Anihilation operator for different magnetic projections (sigma)
ComplexMatrix Annihilation(int sigma){
//Sigma {-1,0,1} -> {0,1,2}
	ComplexMatrix O(m,m);
	for(int i=0;i<m/4;i++){
		for(int j=0;j<m/4;j++){
			
			for(int alpha=0;alpha<m/4;alpha++){
				for(int beta=0;beta<m/4;beta++){
					int x = i*(m/4)+alpha;
					int y = j*(m/4)+beta;
					if(i==j){
						if(alpha==0 && beta>0 ){
							O(x,y)=delta(beta-1,sigma);
							}
						}
					if(i==0 && j>0){
						if(alpha==beta){
							O(x,y)=delta(j-1,sigma);
							}
						}
					}
				}
			}	
		}
	return O;			
}

//Number operator for each magnetic projection
ComplexMatrix small_number_Op(int sigma){
//Sigma {-1,0,1} -> {0,1,2}
	ComplexMatrix B(m,m), BT(m,m),N_matrix(m,m);
	B=Creation(sigma);
	BT= Annihilation(sigma);
	
	N_matrix=prod(B,BT);
	
	return N_matrix;
}

//Complete number operator (calculated with creation and anihilation operators)
ComplexMatrix Number_Op1(){
	ComplexMatrix N_mtrx(m,m);
	ComplexMatrix n_mtrx(m,m);
	
	for(int i=0;i<m/4;i++){
		n_mtrx = small_number_Op(i);
		N_mtrx+= n_mtrx;
	
	}
	
	return N_mtrx;
}

//Complete number operator (manually specified)
ComplexMatrix Number_Op2(){
	ComplexMatrix N_mtrx(m,m);
	for(int l=0;l<3;l++){
	for(int i=5;i<8;i++){
		N_mtrx(i+4*l,i+4*l)= 2;
		}
			}
	for(int i=1;i<4;i++){
	N_mtrx(1*i,1*i)= 1.;
	N_mtrx(4*i,4*i)= 1.;
	}
	return N_mtrx;
}

//General number operator
ComplexMatrix Number_Op(int sig){
//Sigma {-1,0,1} -> {0,1,2} (sig=3 for total number of particles)
	ComplexMatrix N_mtrx(m,m);
	if(sig==3){
		N_mtrx=Number_Op2();
		}
	else{
	N_mtrx=small_number_Op(sig);
		}
	return N_mtrx;
}


//----------------- Functions ------------------//


// Forcing the matrices for S=1 and S=L to be vectors.
 void Correct(ComplexMatrix& A,bool Tail){
  ComplexMatrix U_(M*m,M);
  if(Tail==0){
    SubComplexMatrix Ak(A,Range(0,m*M),Range(0,M));
    for(int Ns=0;Ns<m;Ns++)
      for(int alpha=0;alpha<M;alpha++)
		U_(Ns*M,alpha)=Ak(Ns*M,alpha);
    Ak=U_;
  }else{
    SubComplexMatrix Ak(A,Range(0,m*M),Range((L-1)*M,L*M));
    for(int Ns=0;Ns<m;Ns++)
      for(int alpha=0;alpha<M;alpha++)
		U_(Ns*M+alpha,0)=Ak(Ns*M+alpha,0);
    Ak=U_;
  }
}


//SVD contraction protocol for orthogonalization
void OrthoSD(ComplexMatrix &A,int& Sx,bool dir){
	
  // Definitions which are the same in all cases
  ComplexMatrix Sigma(M,M); 
  RealVector D(M);
  int alpha;
  int Ns;

//  std::cout<<"Case dir="<<dir<<" and Sx="<<Sx<<".\n";// OUTPUT
	
  if(dir){
    // dir=1: Right direction. Sx<L-1 assumed.
    
    
    
    if(Sx==0){
      SubComplexMatrix Ak(A,Range(0,m*M),Range(Sx*M,(Sx+1)*M));
      ComplexMatrix Ak_(m,M);
      ComplexMatrix U(m,M);
      ComplexMatrix V(M,M);

      	// Write first rows of A[1]^{n_1} into Ak_
      for(Ns=0;Ns<m;Ns++)
		for(alpha=0;alpha<M;alpha++)
		  Ak_(Ns,alpha)=Ak(Ns*M,alpha);

      	// Taking the SVD for A(Sx)
      svd(Ak_,U,D,V);

  //    std::cout<<"D("<<Sx<<")="<<D<<"\n"; //OUTPUT

      	// Write U into new A[1] assuming all other entries have been '0' before
      for(Ns=0;Ns<m;Ns++)
		for(alpha=0;alpha<M;alpha++)
		  Ak(Ns*M,alpha)=U(Ns,alpha);    
	  
	// Put the inperfection on Sx+1  
      for(alpha=0;alpha<M;alpha++) 
		Sigma(alpha,alpha)=D(alpha);
	  Sigma=prod(Sigma,V);
	  
      ComplexMatrix ANs_(M,M);     
      for(Ns=0;Ns<m;Ns++){
		SubComplexMatrix ANs(A,Range(Ns*M,(Ns+1)*M),Range((Sx+1)*M,(Sx+2)*M));
		ANs_=ANs;
		noalias(ANs)=prod(Sigma,ANs_);
      }

	//  std::cout<<"Finito case Sx="<<Sx<<".\n";// OUTPUT
	  
    }	// end of case (dir && (Sx==0))
	else{     
      SubComplexMatrix Ak(A,Range(0,m*M),Range(Sx*M,(Sx+1)*M));
      ComplexMatrix Ak_(m*M,M);
      ComplexMatrix U(m*M,M);
      ComplexMatrix V(M,M);

 	// Taking the SVD for A(Sx)  
      Ak_=Ak;
      svd(Ak_,U,D,V);
	  
//	  std::cout<<"D("<<Sx<<")="<<D<<"\n"; // OUTPUT
	  
 	Ak=U;
        ComplexMatrix ANs_(M,M);  
	 // Put the inperfection on Sx+1  
      for(alpha=0;alpha<M;alpha++){
		Sigma(alpha,alpha)=D(alpha);   
      		Sigma=prod(Sigma,V);
	  	ComplexMatrix ANs_(M,M);
	}    
	  	 
      for(Ns=0;Ns<m;Ns++){
		SubComplexMatrix ANs(A,Range(Ns*M,(Ns+1)*M),Range((Sx+1)*M,(Sx+2)*M));
		ANs_=ANs;
		noalias(ANs)=prod(Sigma,ANs_);
      }

//	  std::cout<<"Finito case Sx="<<Sx<<".\n";// OUTPUT
    }// end of case (dir && (Sx>0))
    
    
    Sx++; // Now Sx has moved by one position to the right
   
    
  }else{
    // dir==0: Moving to the left, Sx>1 assumed.
    if(Sx==L-1){	
      SubComplexMatrix Ak(A,Range(0,m*M),Range(Sx*M,(Sx+1)*M));
      ComplexMatrix Ak_(M,m);
      ComplexMatrix U(M,M);
      ComplexMatrix V(M,m);

//	  std::cout<<"Prepare Ak_.\n";// OUTPUT
	  
	  // copy first columns of A[L-1]^n_{L-1} into Ak
      for(Ns=0;Ns<m;Ns++)
		for(alpha=0;alpha<M;alpha++)
		  Ak_(alpha,Ns)=Ak(Ns*M+alpha,0);

//	  std::cout<<"Done.\nPerform SVD.\n";// OUTPUT
	  
      // SVD for the matrix A(Sx)
      svd(Ak_,U,D,V);

//	  std::cout<<"Done.\nD("<<Sx<<")="<<D<<"\n";// OUTPUT
	  
//	  std::cout<<"Write new Ak.\n";// OUTPUT
	  
      for(Ns=0;Ns<m;Ns++)
		for(alpha=0;alpha<M;alpha++)
	  		Ak(Ns*M+alpha,0)=V(alpha,Ns);      
	  
      // Put the imperfection on the left neighbor
      for(alpha=0;alpha<M;alpha++)
		Sigma(alpha,alpha)=D(alpha);
	  Sigma=prod(U,Sigma);

	  ComplexMatrix Ak_l_(m*M,M);
	  SubComplexMatrix Ak_l(A,Range(0,m*M),Range((Sx-1)*M,Sx*M));
	  Ak_l_=Ak_l;
	  noalias(Ak_l)=prod(Ak_l_,Sigma);

//	  std::cout<<"Finito case Sx="<<Sx<<".\n";// OUTPUT
    }// end of case (dir==0 && (Sx==L-1))
	else{
      ComplexMatrix Ak_(M,m*M);
      ComplexMatrix U(M,M);
      ComplexMatrix V(M,m*M);
	 
	  //set Ak_ from A
      for(Ns=0;Ns<m;Ns++){
		SubComplexMatrix Ak_Ns(A,Range(Ns*M,(Ns+1)*M),Range(Sx*M,(Sx+1)*M));
		SubComplexMatrix Ak_Ns_(Ak_,Range(0,M),Range(Ns*M,(Ns+1)*M));
		Ak_Ns_=Ak_Ns;
      }
	  
      // SVD for the matrix A(Sx)
      svd(Ak_,U,D,V);

//	  std::cout<<"D("<<Sx<<")="<<D<<"\n"; // OUTPUT
	 
	  // set new Ak from V
      for(Ns=0;Ns<m;Ns++){
		SubComplexMatrix Ak_Ns(A,Range(Ns*M,(Ns+1)*M),Range(Sx*M,(Sx+1)*M));
		SubComplexMatrix V_Ns(V,Range(0,M),Range(Ns*M,(Ns+1)*M));
		Ak_Ns=V_Ns;
      }
	  
      // Put the imperfection on the left neighbor
      for(alpha=0;alpha<M;alpha++){
		Sigma(alpha,alpha)=D(alpha);
	  	Sigma=prod(U,Sigma);

		ComplexMatrix Ak_l_(m*M,M);
		SubComplexMatrix Ak_l(A,Range(0,m*M),Range((Sx-1)*M,Sx*M));
		Ak_l_=Ak_l;
		noalias(Ak_l)=prod(Ak_l_,Sigma);
	}
//	  std::cout<<"Finito case Sx="<<Sx<<".\n";// OUTPUT
    }// end of case (dir==0 && (Sx<L-1))
    Sx--; // Now the imperfection is in Sx--
  }// end of case (dir==0)
}//---End OrthoSD---//


// This perform the SVD on the left and right hand sides of the site k.
// So, the imperfection remains on the site k.
void Orthogonalize(ComplexMatrix &A,int k){
  int Sx=0;
  while(Sx<k){
    OrthoSD(A,Sx,1);
  }
  Sx=L-1;
  while(Sx>k){
    OrthoSD(A,Sx,0);    
  }
}

//Product of an operator in matrix representation with the matrix A
ComplexMatrix opten(ComplexMatrix &O, ComplexMatrix &A){
  int Od=O.size1() , Ar=A.size1(), Ac=A.size2();
  //Op. dimension    #rows in A    #columns in A
  ComplexMatrix OA(Ar,Ac);
  for(int n=0;n<Od;n++){
      SubComplexMatrix AOn(OA,Range(n*M,(n+1)*M),Range(0,Ac));
    for(int n_=0;n_<Od;n_++){
      SubComplexMatrix An_(A,Range(n_*M,(n_+1)*M),Range(0,Ac));
      AOn+=O(n,n_)*An_;
    }
  }
  return OA;
}

//Function to compute single site operator product with tensor
void Single_site_opten(ComplexMatrix& O,ComplexMatrix& A, int k){
		//A(m*M,L*M) matrix orthonormalized	O operator (m,m)	k site
	int A_r = A.size1();
	SubComplexMatrix Ak(A,Range(0, A_r),Range(k*M, (k+1)*M));
	ComplexMatrix Ak_(A_r, M);
	Ak_=Ak;
	Ak=opten(O,Ak_);
}

//Trace over left and right (d) transfer matrix
ComplexMatrix Trmu(ComplexMatrix &A,bool d){
  int alpha,alpha_,beta;
  ComplexMatrix a(M,M);
  if(d){
    for(alpha=0;alpha<M;alpha++)
      for(alpha_=0;alpha_<M;alpha_++)
	for(beta=0;beta<M;beta++)
	  a(alpha,alpha_)+=A(alpha*M+alpha_,beta*M+beta);
  }else{
  for(alpha=0;alpha<M;alpha++)
    for(alpha_=0;alpha_<M;alpha_++)
      for(beta=0;beta<M;beta++)
	a(alpha,alpha_)+=A(beta*M+beta,alpha*M+alpha_);
    }
  return a;
}

//Transfer matrix for sites without operator
  ComplexMatrix TrnsfrMtrx(ComplexMatrix &A, int k){
  int alpha,beta,Ns,alpha_,beta_;
  SubComplexMatrix Ak_(A,Range(0,m*M),Range(k*M,(k+1)*M));
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
  
//Transfer matrix for sites with operator
ComplexMatrix TrnsfrMtrx_(ComplexMatrix &A,ComplexMatrix &B, int k){
  int alpha,beta,Ns,alpha_,beta_;
  SubComplexMatrix Ak_(A,Range(0,m*M),Range(k*M,(k+1)*M));
  SubComplexMatrix Bk_(B,Range(0,m*M),Range(k*M,(k+1)*M));
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

//State norm calculation
Complex Norm_(ComplexMatrix& A,int site_information){
		Complex norm_;
		int Ar=A.size1(), Ac=A.size2();
		ComplexMatrix A0(Ar,Ac); // m*M, M*L
		
		ComplexMatrix O(m,m);
		

		Orthogonalize(A,site_information);
		
		A0 = A; //a copy of the ket->future bra
		
		//Usando la función opten para calcular el operado por matriz de un sitio
                O=ident_op();
		
		Single_site_opten(O,A,site_information);
		
		
		ComplexVector E(M*M),Ef(M*M);  
		ComplexMatrix EE_1(M,M),EE_2(M*M,M*M),EE_3(M*M,M*M);


		EE_2=TrnsfrMtrx_(A0,A,0);
		EE_1=Trmu(EE_2,0);

		for(int alpha=0;alpha<M;alpha++)
			for(int alpha_=0;alpha_<M;alpha_++)
			  E(alpha*M+alpha_)=EE_1(alpha,alpha_);

		for(int ss=1;ss<L;ss++){
			if(ss<L-1){
			  EE_3= TrnsfrMtrx_(A0,A,ss);
			  E=prod(trans(EE_3),E);
			}
		
		else{	
			EE_2=TrnsfrMtrx_(A0,A,ss);
			EE_1=Trmu(EE_2,1);
			for(int alpha=0;alpha<M;alpha++)
			  for(int alpha_=0;alpha_<M;alpha_++)
			    Ef(alpha*M+alpha_)=EE_1(alpha,alpha_);
			}
		}
		
		norm_=inner_prod(E,Ef); 
	
	return norm_;
}

//Calculates the expected value of an operator on a given site 
Real Single_site_op(ComplexMatrix& A,ComplexMatrix&Op, int site){
		Real measure;
		int Ar=A.size1(), Ac=A.size2();

		//cout << "Ar= "<< A.size1() << endl;
		//cout << "Ac= "<< A.size2() << endl;		
		
		ComplexMatrix Abra(Ar,Ac),Aket(Ar,Ac); // m*M, M*L
		Abra = A; //a copy of the ket->future bra
		Aket = A; //a copy of the ket
		
		ComplexMatrix O(m,m);
		
		O=Op;
		
		Single_site_opten(O,Aket,site);
		
		ComplexVector E(M*M),Ef(M*M);  
		ComplexMatrix EE_1(M,M),EE_2(M*M,M*M),EE_3(M*M,M*M);

		EE_2=TrnsfrMtrx_(Abra,Aket,0);
		EE_1=Trmu(EE_2,0);

		for(int alpha=0;alpha<M;alpha++)
			for(int alpha_=0;alpha_<M;alpha_++)
			  E(alpha*M+alpha_)=EE_1(alpha,alpha_);

		for(int ss=1;ss<L;ss++){
			if(ss<L-1){
			  EE_3= TrnsfrMtrx_(Abra,Aket,ss);
			  E=prod(trans(EE_3),E);
			}
		else{	
			EE_2=TrnsfrMtrx_(Abra,Aket,ss);
			EE_1=Trmu(EE_2,1);
			for(int alpha=0;alpha<M;alpha++)
			  for(int alpha_=0;alpha_<M;alpha_++)
			    Ef(alpha*M+alpha_)=EE_1(alpha,alpha_);
			}
		}
		measure=real(inner_prod(E,Ef)); 
	return measure;
}
// This routine makes a othogonalizacion puting the imperfection in the site k
// and then it normalizes the wave fuction represented by the MPS.
void Orthonormalize1(ComplexMatrix &A,int k){
  Real N=0;
  
  Orthogonalize(A,k);
  
  for(int ns=0;ns<m;ns++)
    for(int alpha=0;alpha<M;alpha++)
      for(int alpha_=0;alpha_<M;alpha_++)
	N+=real(conj(A(ns*M+alpha,k*M+alpha_))*A(ns*M+alpha,k*M+alpha_));
  //std::cout<<"N="<<N<<"\n";
  for(int ns=0;ns<m;ns++)
    for(int alpha=0;alpha<M;alpha++)
      for(int alpha_=0;alpha_<M;alpha_++)
	A(ns*M+alpha,k*M+alpha_)/=sqrt(N);
}

//Orthonormalization using norm calculated using transfer matrices 
void Orthonormalize2(ComplexMatrix& A, int k){
	ComplexMatrix Op(m,m);
	Op=ident_op();
	Real norm=Single_site_op(A,Op,k);

	for(int ns=0;ns<m;ns++){
		for(int alpha=0;alpha<M;alpha++){
			for(int alpha_=0;alpha_<M;alpha_++){
				A(ns*M+alpha,k*M+alpha_)/=sqrt(norm);
			}
		}
	}
}


