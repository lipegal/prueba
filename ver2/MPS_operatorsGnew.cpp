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
		if(S==1) A(5*M,S*M)=1.;
		else if(S==3) A(10*M,S*M)=1.; 
		else if(S==2) A(15*M,S*M)=1.;
		else if(S==4) A(8*M,S*M)=1.;
		else if(S==5) A(13*M,S*M)=1.;
		else if(S==6) A(7*M,S*M)=1.;
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

int main(){


FILE *dskw;
char archive[300];
snprintf(archive, sizeof(archive), "output.txt");
dskw=fopen(archive,"w+");


ComplexMatrix A(M*m,M*L); 
Complex norm;


A_mat(A);
//RACalc(A);
//ACalc(A);
//cout<<A<<endl; 

int site_information = 2;
int k= site_information;


cout<<"Comienza"<<endl;

norm=Norm_(A,k);

cout<<"norma antes: "<<real(norm)<<" "<<imag(norm)<<endl;

Orthonormalize2(A,k);

cout<<"norma después: "<<real(norm)<<" "<<imag(norm)<<endl;

Real Particles=0.;
Real Total_Particles=0.;

for(int l=0; l<L;l++){
	Particles=Number_Particles(A,l,3);
	cout<<"Particles site "<<l<<": "<<Particles<<endl;
	Total_Particles+=Particles;
}

cout<<"Total particles "<<Total_Particles<<endl<<endl;
cout<<endl;

for(int m=0;m<3;m++){
Real Particles=0.;
Real Total_Particles=0.;
cout<<"Magnetic projection= "<<(m-1)<<endl;
	for(int l=0; l<L;l++){
		Particles=Number_Particles(A,l,m);
		cout<<"Particles site "<<l<<": "<<Particles<<endl;
		Total_Particles+=Particles;
}
cout<<"Total particles="<<Total_Particles<<endl<<endl;
cout<<"--------------------------------------------"<<endl;
}



//fprintf(dskw,"---------------------------------\n");

fflush(dskw);
return 0;
}



