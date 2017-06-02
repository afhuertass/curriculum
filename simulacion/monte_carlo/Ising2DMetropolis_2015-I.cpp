#include<iostream>
#include<cmath>
#include "Random64.h"
using namespace std;

const int L=12;
const int L2=L*L;

class SpinSystem{
private:
  int s[L][L], E,M;
public:
  void InicieArriba(void);
  double GetE(void){return (double) E;};
  double GetM(void){return fabs((double)M);};
  void UnPasoDeMetropolis(Crandom &ran2,double Beta);
};
void SpinSystem::InicieArriba(void){
  for(int i=0;i<L;i++)
    for(int j=0;j<L;j++)
      s[i][j]=1;
  E=-2*L2; M=L2;
}
void SpinSystem::UnPasoDeMetropolis(Crandom &ran2,double Beta){
  int n,i,j,dE;
  //Escoger un espin al azar;
  n=(int) (L2*ran2.r()); i=n/L; j=n%L;
  //Calcular dE;
  dE=2*s[i][j]*(s[i][(j+1)%L]+s[i][(j-1+L)%L]+s[(i+1)%L][j]+s[(i-1+L)%L][j]);
  //Aceptar segun la regla de Metropolis A(x'|x)
  if(dE<=0)
    {s[i][j]*=-1; E+=dE; M+=2*s[i][j];} //Acepto voltear el espin
  else if(ran2.r()<exp(-Beta*dE))
    {s[i][j]*=-1; E+=dE; M+=2*s[i][j];} //Acepto voltear el espin
}

int main(void){

  const int EquilibriumSteps=(int)(200*pow(L/8.0,2.125)); //MCSS
  const int CorrSteps=(int)(50*pow(L/8.0,2.125)); //MCSS
  const int NMuestras=10000;
  const double k=1;

  SpinSystem Ising;
  Crandom ran2(10);
  int mcs,s_eq,s_corr,n; double kT;
  double SumM,SumM2,SumM4,SumE,SumE2,E,M;
  double Mprom,M2prom,M4prom,Eprom,E2prom,Cv,Xs,Ubinder;

  for(kT=0.4;kT<4.0;kT+=0.2){
    //INICIALIZAR;
    Ising.InicieArriba();
    //LLEGAR AL EQUILIBRIO;
    for(s_eq=0;s_eq<EquilibriumSteps;s_eq++)
      for(mcs=0;mcs<L2;mcs++) //un paso de Monte Carlo por sitio (MCSS)
	Ising.UnPasoDeMetropolis(ran2,1.0/kT);
    
    //TOMAR NMUESTRAS
    SumM=SumM2=SumM4=SumE=SumE2=0; //Arranque acumuladores en 0
    for(n=0;n<NMuestras;n++){
      //TOMAR LOS DATOS;
      E=Ising.GetE();  M=Ising.GetM();
      SumM+=M; SumM2+=M*M; SumM4+=M*M*M*M; SumE+=E; SumE2+=E*E; 
      //AVANZAR A LA SIGUIENTE MUESTRA;
      for(s_corr=0;s_corr<CorrSteps;s_corr++)
	for(mcs=0;mcs<L2;mcs++) //un paso de Monte Carlo por sitio (MCSS)
	  Ising.UnPasoDeMetropolis(ran2,1.0/kT);
    }
    
    //POST-PROCESAR
    Mprom=SumM/NMuestras; M2prom=SumM2/NMuestras; M4prom=SumM4/NMuestras;
    Eprom=SumE/NMuestras; E2prom=SumE2/NMuestras;
    
    Cv=k/(kT*kT)*(E2prom-Eprom*Eprom);
    Xs=1.0/kT*(M2prom-Mprom*Mprom);
    Ubinder=1-1.0/3*(M4prom/(M2prom*M2prom));
    
    //IMPRIMIR RESULTADOS
    cout<<kT<<" "<<Eprom<<" "<<Mprom<<" "<<Cv<<" "<<Xs<<" "<<Ubinder<<endl;
  }
  return 0;
}
