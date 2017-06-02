// Difusion 1D 
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=128;
const int Ly=64;

const double Tau=2.85;
const double RHO0=1.0, UX0=0.03, UY0=0;

const int ixc=50, iyc=32, R=12;

enum TipoCelda{aire,obstaculo,ventilador};

class LatticeBoltzmann{
private:
  double f[Lx][Ly][9],fnew[Lx][Ly][9];//f[ix][iy][i]
  int V[2][9]; // V[x=0,y=1][i]
  double w[9]; // w[i]
  TipoCelda Celda[Lx][Ly];
public:
  LatticeBoltzmann(void);
  void Inicie(void);
  void ConstruyaLaGeometria(void);
  double rho(int ix,int iy);
  double Ux(int ix,int iy);
  double Uy(int ix,int iy);
  double feq(double rho0,double Ux0,double Uy0,int i);
  void Colisione(void);
  void Adveccione(void);
  void Imprimase(char * NombreArchivo);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //Cargar los pesos
  w[0]=4.0/9;    w[1]=w[2]=w[3]=w[4]=1.0/9;    w[5]=w[6]=w[7]=w[8]=1.0/36;
  //Cargar los vectores
  V[0][0]=0;  
  V[1][0]=0;
  
  V[0][1]=1;  V[0][2]=0;  V[0][3]=-1; V[0][4]=0;  
  V[1][1]=0;  V[1][2]=1;  V[1][3]=0;  V[1][4]=-1;

  V[0][5]=1;  V[0][6]=-1; V[0][7]=-1; V[0][8]=1;  
  V[1][5]=1;  V[1][6]=1;  V[1][7]=-1; V[1][8]=-1;
}
void LatticeBoltzmann::Inicie(void){
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<9;i++)
	fnew[ix][iy][i]=f[ix][iy][i]=feq(RHO0,UX0,UY0,i);
}
void LatticeBoltzmann::ConstruyaLaGeometria(void){
  int ix,iy;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      if(ix==0)
	Celda[ix][iy]=ventilador;
      else if((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)<=R*R)
	Celda[ix][iy]=obstaculo;
      else
	Celda[ix][iy]=aire;
}
double LatticeBoltzmann::rho(int ix,int iy){
    int i; double suma;
    for(suma=0,i=0;i<9;i++)
      suma+=f[ix][iy][i];
    return suma;
}
double LatticeBoltzmann::Ux(int ix,int iy){
  if(Celda[ix][iy]==ventilador)
    return UX0;
  else if(Celda[ix][iy]==obstaculo)
    return 0;
  else{
    int i; double suma;
    for(suma=0,i=0;i<9;i++)
      suma+=V[0][i]*f[ix][iy][i];
    return suma/rho(ix,iy);
  }
}
double LatticeBoltzmann::Uy(int ix,int iy){
  if(Celda[ix][iy]==ventilador)
    return UY0;
  else if(Celda[ix][iy]==obstaculo)
    return 0;
  else{
    int i; double suma;
    for(suma=0,i=0;i<9;i++)
      suma+=V[1][i]*f[ix][iy][i];
    return suma/rho(ix,iy);
  }
}
double LatticeBoltzmann::feq(double rho0,double Ux0,double Uy0,int i){
  double U2=Ux0*Ux0+Uy0*Uy0; double UpVi=Ux0*V[0][i]+Uy0*V[1][i];
  return rho0*w[i]*(4.5*UpVi*UpVi+3*UpVi+1-1.5*U2);
}
void LatticeBoltzmann::Colisione(void){
  int ix,iy,i; double rho0,Ux0,Uy0;
  for(ix=0;ix<Lx;ix++)//Para cada celda
    for(iy=0;iy<Ly;iy++){
      //Calcular las cantidades macroscopicas
      rho0=rho(ix,iy);   Ux0=Ux(ix,iy);   Uy0=Uy(ix,iy);
      for(i=0;i<9;i++)//Para cada direccion en la celda
	//Hago la evolución de la Ec. de Boltzmann BGK
	fnew[ix][iy][i]=f[ix][iy][i]-1/Tau*(f[ix][iy][i]-feq(rho0,Ux0,Uy0,i));
    }
}
void LatticeBoltzmann::Adveccione(void){
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<9;i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
}
void LatticeBoltzmann::Imprimase(char * NombreArchivo){
  ofstream MiArchivo(NombreArchivo); int ix,iy;
  for(ix=0;ix<Lx;ix+=4)
    for(iy=0;iy<Ly;iy+=4)
      MiArchivo<<ix<<" "<<iy<<" "<<20*Ux(ix,iy)<<" "<<20*Uy(ix,iy)<<endl;
  MiArchivo.close();
}

const int tmax=200;

int main(){
  LatticeBoltzmann Aire;
  int t;
  
  Aire.ConstruyaLaGeometria();
  Aire.Inicie();
  
  for(t=0;t<tmax;t++){
    Aire.Adveccione();
    Aire.Colisione();
  }
  Aire.Imprimase("Aire.dat");
  
  return 0;
}
