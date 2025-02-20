
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iomanip>
#include <iostream>
#include <lapack.h>
#include <matrixtypes.h>

using std::cout;
using std::cin;
using std::endl;
using std::flush;

using namespace ula;

#define PI 3.14159265358979323846

int m=16,M=20,L=10;

#include <MPSlib.h>

void A_mat(ComplexMatrix& A){
	for(int S=0;S<L;S++){
		if(S==3){
			A(11*M,S*M)=sqrt(2);
			A(14*M,S*M)=sqrt(2);
		}
		//if(S==3) A(3*M,S*M)=1.; 
		//if(S==1) A(3*M,S*M)=1.; 
		//if(S==5){ A(3*M,S*M)=0.5; A(5*M,S*M)=0.5;} 
		else A(0,S*M)=1.; //No particles on other sites 
	}
	Correct(A,0); 
	Correct(A,1);
}

void ACalc(ComplexMatrix& A){
  int S;
  for(S=0;S<L;S++)
    //    A(0,S*M)=1.;
    { 
      A(0,S*M)  =sqrt(1.-2.*Gauss(S,10,1.5));
      A(M,S*M)  =sqrt(2.*Gauss(S,10,1.5))*exp(I*(double)S*PI/2.);
      A(2*M,S*M)=sqrt(2.*Gauss(S,10,1.5))*exp(-I*(double)S*3.*PI/2.);
    }
    /*  if(S==20)
    A(2*M,S*M)=1.;
  else
  A(0,S*M)=1.; */
  Correct(A,0);
  Correct(A,1); 
}



ComplexMatrix ident_op(){
  ComplexMatrix Id(m,m);
  for(int n=0;n<m;n++)
    Id(n,n)=1.;
  return Id;  
}

//Applies an operator on a single site using opten
void Single_Site_Operator_Procedure(ComplexMatrix& O,ComplexMatrix& A, int k){
			//A(m*M,L*M) matrix orthonormalized	O operator (m,m)	k site
	int A_r = A.size1();
	SubComplexMatrix Ak(A,Range(0,A_r),Range(k*M,(k+1)*M));
	ComplexMatrix Ak_(A_r,M);

		cout << "Ak_r= "<< Ak.size1() << endl;
		cout << "Ak_c= "<< Ak.size2() << endl;
		
		cout << "Ak__r= "<< Ak_.size1() << endl;
		cout << "Ak__c= "<< Ak_.size2() << endl;
	Ak_=Ak;
	Ak=opten(O,Ak_);
	
}


Complex Norm_(ComplexMatrix& A,int k){
		Complex norm_;
		int A_r=A.size1(), A_c=A.size2();

		ComplexMatrix A0(A_r,A_c); // m*M, M*L
		ComplexMatrix O(m,m);

		Orthogonalize(A,k);
		
		A0 = A; //a copy of the ket->future bra
		
		//Usando la función opten para calcular el operador por matriz de un sitio
        O=ident_op();
		//Number_Op();
		
		Single_Site_Operator_Procedure(O,A,k);

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
			//cout<<"Reached, L="<<ss<<endl;

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


Real Number_Particles(ComplexMatrix& A,int k){
		Real number;
		int A_r=A.size1(), A_c=A.size2();

		cout << "Ar= "<< A.size1() << endl;
		cout << "Ac= "<< A.size2() << endl;		
		
		ComplexMatrix Abra(A_r,A_c),A_ket(A_r,A_c); // m*M, M*L

		Orthogonalize(A,k);
		Abra = A; //a copy of the ket->future bra
		A_ket = A; //a copy of the ket
		
		ComplexMatrix O(m,m);

		O=Number_Op();
		
		Single_Site_Operator_Procedure(O,A_ket,k);
		
		ComplexVector E(M*M),Ef(M*M);  
		ComplexMatrix EE_1(M,M),EE_2(M*M,M*M),EE_3(M*M,M*M);

		EE_2=TrnsfrMtrx_(Abra,A_ket,0);
		EE_1=Trmu(EE_2,0);

		for(int alpha=0;alpha<M;alpha++)
			for(int alpha_=0;alpha_<M;alpha_++)
			  E(alpha*M+alpha_)=EE_1(alpha,alpha_);

		for(int ss=1;ss<L;ss++){
			if(ss<L-1){
			  EE_3= TrnsfrMtrx_(Abra,A_ket,ss);
			  E=prod(trans(EE_3),E);
			}
		else{	
		//cout<<"Reached L= "<<ss<<endl;
			EE_2=TrnsfrMtrx_(Abra,A_ket,ss);
			EE_1=Trmu(EE_2,1);
			for(int alpha=0;alpha<M;alpha++)
			  for(int alpha_=0;alpha_<M;alpha_++)
			    Ef(alpha*M+alpha_)=EE_1(alpha,alpha_);
			}
		}
		
		number=real(inner_prod(E,Ef)); 
	
	return number;
}

int main(){


FILE *dskw;
char archive[300];
snprintf(archive, sizeof(archive), "output.txt");
dskw=fopen(archive,"w+");



ComplexMatrix A(M*m,M*L), A0(M*m,M*L); 
Complex norm;


A_mat(A);
//RACalc(A);
//ACalc(A);
//cout<<A<<endl; 

int site_information = 3;
int k= site_information;


norm=Norm_(A,k);
cout<<"norma before: "<<real(norm)<<" "<<imag(norm)<<endl;

fprintf(dskw,"(Norm) %d %f %f\n",M,real(norm),imag(norm));

for(int ns=0;ns<m;ns++)
    for(int alpha=0;alpha<M;alpha++)
      for(int alpha_=0;alpha_<M;alpha_++)
	A(ns*M+alpha,k*M+alpha_)/=sqrt(real(norm));
	
norm=Norm_(A,k); //Normalized

fprintf(dskw,"(Norm) %d %f %f\n",M,real(norm),imag(norm));

cout<<"norm after: "<<real(norm)<<" "<<imag(norm)<<endl;

/*
Real Nn=0.;

 for(int ns=0;ns<m;ns++)
    for(int alpha=0;alpha<M;alpha++)
      for(int alpha_=0;alpha_<M;alpha_++)
	A(ns*M+alpha,k*M+alpha_)/=sqrt(Nn);
*/

Complex Particles=0.;

Particles=Number_Particles(A,k);

cout<<"Particle number: "<<Particles<<endl;

fprintf(dskw,"(Number of particles) %d %f \n",M,real(Particles));
fprintf(dskw,"---------------------------------\n");
//cout<<Number_Op();
fflush(dskw);
return 0;
}



