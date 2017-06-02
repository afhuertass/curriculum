// Programa que simula la cuna de newton creando una 
// animación en gnuplot
#include <iostream>
#include <cmath>
using namespace std;

const double g=980;
const double Deltat=  1e-6;
const double chi=0.193183325037836;

const int N=5;
const double K=1e8;

class Cuerpo{
private:
  double theta,omega,tau; double R,L,I,m,x0;
public:
  void Inicie(double theta0,double omega0,double m0,double R0,double L0,double x00);
  void CalculeTorqueBase(void);
  void CalculeTorqueQueMeHace(Cuerpo ElOtroPendulo);
  void Mueva_r1(double dt);
  void Mueva_V(double dt);
  void Mueva_r2(double dt);
  void Dibujese(void);
  double Getx(void){return x0+L*sin(theta);};//Funcion inline
  double Gety(void){return -L*cos(theta);};//Funcion inline
  double GetTau(void){return tau;};//Funcion inline
};
void Cuerpo::Inicie(double theta0,double omega0,double m0,double R0,double L0,double x00){
  theta=theta0; omega=omega0; m=m0; R=R0; L=L0; x0=x00; I=m*L*L;
}
void Cuerpo::CalculeTorqueBase(void){
  tau=-L*m*g*sin(theta);
}
void Cuerpo::CalculeTorqueQueMeHace(Cuerpo ElOtroPendulo){
  double dx=Getx()-ElOtroPendulo.Getx();
  double h=(R+ElOtroPendulo.R)-fabs(dx);
  if(h>0) tau+=dx/fabs(dx)*K*pow(h,1.5);
}
void Cuerpo::Mueva_r1(double dt){
  theta+=omega*chi*dt;
}
void Cuerpo::Mueva_V(double dt){
  omega+=tau*dt/(2*I);
}
void Cuerpo::Mueva_r2(double dt){
  theta+=omega*(1-2*chi)*dt;
}
void Cuerpo::Dibujese(void){
  //cout<<", "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+"<<R<<"*sin(t)";
 cout<<", "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+"<<R<<"*sin(t)";
}

//----------- Funciones Globales --------
void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl; 
  //cout<<"set output 'UnPendulo.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-10:18]"<<endl;
  cout<<"set yrange [-14:0]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}

int main(){
  Cuerpo Pendulo[N];
  double t,tdibujo;
  double R0=1,M0=20,L0=10;
  double T=10; int i;
  
  T=2*M_PI*sqrt(L0/g);
  
    InicieAnimacion();
  //------------(theta0     ,omega0,m0,R0,L0,x00)
  Pendulo[0].Inicie(-15*M_PI/180,0     ,M0,R0,L0,0  );
  for(i=1;i<N;i++) Pendulo[i].Inicie(0,0     ,M0,R0,L0,2*R0*i);
  for(t=tdibujo=0;t<T;t+=Deltat,tdibujo+=Deltat){
    
    
      
      if(tdibujo > T/100) {
	InicieCuadro();
	for(i=0;i<N;i++) Pendulo[i].Dibujese();
	TermineCuadro();
	tdibujo=0;
      } 
    
    //cout<<t<<" "<<Pendulo[1].GetTau()<<" "<<endl;

    for(i=0;i<N;i++) Pendulo[i].Mueva_r1(Deltat);

    //CalculeTorques
    for(i=0;i<N;i++)   Pendulo[i].CalculeTorqueBase();
    for(i=0;i<N-1;i++) Pendulo[i].CalculeTorqueQueMeHace(Pendulo[i+1]);
    for(i=1;i<N;i++)   Pendulo[i].CalculeTorqueQueMeHace(Pendulo[i-1]);

    for(i=0;i<N;i++){Pendulo[i].Mueva_V(Deltat); Pendulo[i].Mueva_r2(Deltat);}

    //CalculeTorques
    for(i=0;i<N;i++)   Pendulo[i].CalculeTorqueBase();
    for(i=0;i<N-1;i++) Pendulo[i].CalculeTorqueQueMeHace(Pendulo[i+1]);
    for(i=1;i<N;i++)   Pendulo[i].CalculeTorqueQueMeHace(Pendulo[i-1]);

    for(i=0;i<N;i++){Pendulo[i].Mueva_V(Deltat); Pendulo[i].Mueva_r1(Deltat);}
    
  }
  return 0;
}
