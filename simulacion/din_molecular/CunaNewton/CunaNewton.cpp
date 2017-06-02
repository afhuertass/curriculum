// Mi primer programa de DINAMICA MOLECULAR :)

#include <iostream>
#include <cmath>
using namespace std;


const double g=980;
const double Deltat=0.0001;  
const double chi=0.1931833250;
const int N=5;

class Cuerpo{
private: 
  double theta, omega, N; double R,L,I,m,x0;
  
public:
  void Inicie(double theta0, double omega0, double m0, double R0, double L0, double x00);
  void CalculaTorqueBase(void);      
  void Mueva_r1(double dt);
  void Mueva_r2(double dt);
  void Mueva_v(double dt);
  void Dibujese(void);
  double Getx(void){return x0+L*sin(theta);}; //Funcion inline
  double Gety(void){return   -L*cos(theta);};  
}; 

void Cuerpo::Inicie(double theta0, double omega0, double m0, double R0, double L0, double x00){
  theta=theta0; omega=omega0; m=m0; R=R0; L=L0; x0=x00; I=L*L*m;
}

void Cuerpo::CalculaTorqueBase(void){
  N=-L*m*g*sin(theta);
}

void Cuerpo::Mueva_r1(double dt){
  theta+= omega*chi*dt;
}
void Cuerpo::Mueva_v(double dt){
  omega += N*dt/(2*I); 
}
void Cuerpo::Mueva_r2(double dt){
  theta+= omega*(1-2*chi)*dt;
}


void Cuerpo::Dibujese(void){
  cout<<", "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+ "<<R<<"*sin(t)"<<endl;
}




//------FUNCIONES GLOBALES----------------------------

void InicieAnimacion(void){
  cout<<"set terminal gif "<<endl;
  cout<<"set output 'UnPendulo.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-18:18]"<<endl;
  cout<<"set yrange[-14:0]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange[0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;
}

void InicieCuadro(void){
  cout<<"plot 0,0 ";//El espacio es muy importante
}

void TermineCuadro(void){
  cout<<endl;
}

//const double tMAX=10;

//--------Funcion Inicio-----------
int main(){
  Cuerpo Pendulo;
  double t,tdibujo;
  double R0=1,M0=20, L0=10;
  double T;
  T=2*M_PI*sqrt(L0/g);

 
  InicieAnimacion();
  //-------   (theta0      ,omega0,m0,R0,L0,x00);
  Pendulo[0].Inicie(-15*M_PI/180 ,0     ,M0,R0,L0, 0)
  for(i=1;i<N;i++) {
    pendulo[i].Inicie(-15*M_PI/180 ,0     ,M0,R0,L0, 0)
    Pendulo.Inicie(15*M_PI/180 ,0     ,M0,R0,L0, 0);
 }
}
  for(t=tdibujo=0;t<3*T;t+=Deltat,tdibujo+=Deltat){
    
    //cout<<Pendulo.Getx()<<" "<<Pendulo.Gety()<<endl;
    
    if(tdibujo>T/100){
    InicieCuadro();
    Pendulo.Dibujese();
    TermineCuadro();
    tdibujo=0;

    }
    for(i=0; i<N; i++) {Pendulo[i].Mueva_r1(Deltat);}
    for(i=0; i<N; i++) {Pendulo[i].CalculaTorqueBase();}
    for(i=0; i<N; i++) {Pendulo[i].Mueva_v(Deltat);   Pendulo.Mueva_r2(Deltat);}
    for(i=0; i<N; i++) {Pendulo[i].CalculaTorqueBase();}
    for(i=0; i<N; i++) {Pendulo[i].Mueva_v(Deltat); Pendulo.Mueva_r1(Deltat);}

   
  }
  

  return 0;
  
}
