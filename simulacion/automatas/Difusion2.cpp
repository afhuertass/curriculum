// Difusion 1D 
#include <iostream>
#include <cmath>
#include "Random64.h"
using namespace std;
// Automata difusion continua
const int Lx=256;
const double p=0.5;

class Automata{
private:
  double **n, **nnew; //int n[Lx][2],nnew[Lx][2];
  int Vx[2];
public:
  Automata(void);
  ~Automata(void);
  void Inicie(double mu,double sigma);
  void Evolucione();
  double Getn(int ix){return n[ix][0]+n[ix][1];};
  void Muestre(void);
  double GetSigma2(void);
};
Automata::Automata(void){
  n=new double* [Lx];  nnew=new double* [Lx];
  for(int ix=0;ix<Lx;ix++)
    {n[ix]=new double [2];  nnew[ix]=new double [2];}
  Vx[0]=1;  Vx[1]=-1;
}
Automata::~Automata(void){
  for(int ix=0;ix<Lx;ix++)
    {delete[] n[ix];  delete[] nnew[ix];}
  delete[] n;  delete[] nnew;
}
void Automata::Inicie(double mu,double sigma){
  int ix,j;
  //Iniciar todos los contenidos en cero
  for(ix=0;ix<Lx;ix++)
    for(j=0;j<2;j++)
      n[ix][j]=0.5/(sigma*sqrt(2*M_PI))*exp(-0.5*pow((ix-mu)/sigma,2) ) ;
  //Poner la mitad que va a la derecha
 
}
void Automata::Evolucione(){
  int ix,j; double** aux;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(j=0;j<2;j++)
      nnew[(ix+Vx[j]+Lx)%Lx][j]=p*n[ix][j]+ (1-p)*n[ix][(j+1)%2];
 
    aux=n; n=nnew; nnew=aux; //Intercambio los apuntadores
   
}
void Automata::Muestre(void){
  for(int j=0;j<2;j++){
    for(int ix=0;ix<Lx;ix++)
      cout<<n[ix][j];
    cout<<endl;
  }
  cout<<endl;
}

 double Automata::GetSigma2(){
  int ix; double n,ixprom,sigma2;
  //Contar n
  n=0; for(ix=0;ix<Lx;ix++)
	 n+= Getn(ix);
  //Calcular ixprom
  ixprom=0; for(ix=0;ix<Lx;ix++)
	      ixprom+=ix*Getn(ix);
  
  ixprom/=n;
  //Calcular sigma2
  sigma2=0;  for(ix=0;ix<Lx;ix++)
		sigma2+=(ix-ixprom)*(ix-ixprom)*Getn(ix);
  sigma2/=n;

  return sigma2;
}

const int tmax=1000;

int main(){
  Automata Difusion;
  Crandom ran2(10);
  int t;
  
  Difusion.Inicie(0.5*Lx,15.0);
  for ( t = 0; t < tmax ; t++){
    cout << t <<" "  << Difusion.GetSigma2() << endl;
    Difusion.Evolucione();
  }

  return 0;
}
