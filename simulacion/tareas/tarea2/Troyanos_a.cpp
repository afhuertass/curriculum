#include<iostream>
#include<cmath>
#include "Vector.h"
/* ESTE PROGRAMA CORRESPONDE AL PUNTO A DE LA SEGUNDA TAREA
// SIMPLEMENTE SE OBTIENEN LAS COORDENADAS DE JUPITER EN EL SISTEMA NO ROTADO.
Y SE IMPRIMIEN PARA SE GUARDADAS EN UN ARCHIVO.
AL graficarlo se obtiene la trayectoria circular alrededor del sol
*/
using namespace std;

const double Deltat=0.01;
const double G=1;
const double chi=0.193183325037836;
const double Um2chi=1-2*chi;

const int N=3;
//-------------------------------------------------------------
class Cuerpo;
class Colisionador;
//-------------------------------------------------------------

class Cuerpo{
private:
  vector3D r,V,F;  double m,R; //Posicion, Vel, Frueza, masa y Radio

public:
  void Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0);
  void BorreFuerza(void);
  void IncrementeFuerza(vector3D F0);
  void Mueva_r1(double dt);
  void Mueva_V(double dt);
  void Mueva_r2(double dt);
  void Dibujese(void);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  vector3D GetR(void){ return r; };
  friend class Colisionador;
};   // <-------- TIENE QUE COLOCARSE ESE ;!!!!!!!!!!

void Cuerpo::Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0){ //No es una funcion suelta sino que pertence a Cuerpo
  r.cargue(x0,y0,0);  V.cargue(Vx0,Vy0,0);
  m=m0, R=R0;
}
void Cuerpo::BorreFuerza(void){
  F.cargue(0,0,0);
}
void Cuerpo::IncrementeFuerza(vector3D F0){
  F+=F0;
}
void Cuerpo::Mueva_r1(double dt){
  r+=V*(chi*dt); 
}
void Cuerpo::Mueva_V(double dt){
  V+=F*(dt/(2*m));
}
void Cuerpo::Mueva_r2(double dt){
  r+=V*(Um2chi*dt);
}
void Cuerpo::Dibujese(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}


//------------------------------------------------------------------
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Planeta);
  void AgregueFuerza(Cuerpo & Planeta1, Cuerpo & Planeta2);
};

void Colisionador::CalculeFuerzas(Cuerpo * Planeta){
  int i,j;
  for(i=0;i<N;i++)
    Planeta[i].BorreFuerza();
  for(i=0;i<N;i++)
    for(j=0;j<i;j++)
      AgregueFuerza(Planeta[i],Planeta[j]);
}
void Colisionador::AgregueFuerza(Cuerpo & Planeta1, Cuerpo & Planeta2){
  vector3D dr=Planeta1.r-Planeta2.r;
  double aux=G*Planeta1.m*Planeta2.m*pow(norma2(dr),-1.5);
  vector3D F2=dr*aux;
  Planeta2.IncrementeFuerza(F2);  Planeta1.IncrementeFuerza(F2*(-1));
}

//--------------------FUNCIONES GLOBALES-----------------------
void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl;
  //cout<<"set output 'MiPlaneta.gif'"<<endl;
  cout<<""<<endl;
  cout<<"unset key"<<endl;              // QUITAR LABEL
  cout<<"set xrange [-120:120]"<<endl;     // RANGO EN X
  cout<<"set yrange [-120:120]"<<endl;  
  cout<<"set size ratio -1"<<endl;      // VENTANA CUADRADA
  cout<<"set parametric"<<endl;         // GRAFICA PARAMETRICA
  cout<<"set trange [0:7]"<<endl;       // RANGO 
  cout<<"set isosamples 12"<<endl;      // EL BALON ESTA HECHO DE 12 LINEAS

}
void InicieCuadro(void){
  cout<<"plot 0,0 ";
}
void TermineCuadro(void){
  cout<<endl;
}



int main(){
  Cuerpo Planeta[N];
  Colisionador Newton;
  double t,tdibujo;
  double r=100, m0=1047, m1=1;  // datos de entrada para calcular la condcicion inicial
  double r0, r1, V0,V1;
  double M,omega,T;
  int i;
  double mt = 0.005;
  double ang = 60*M_PI/180.0;
  double ang2 = 30*M_PI/180.0;
  // parametros del movimiento y condiciones iniciales.
  
  M=m0+m1;   omega=sqrt(G*M*pow(r,-3));  T=2*M_PI/omega;  
  r1=m0*r/M;  r0=r1-r;   V1=omega*r1;   V0=omega*r0;

  Planeta[0].Inicie(r0, 0 ,0 ,V0 ,m0 , r/10); // sol
  Planeta[1].Inicie(r1, 0 ,0 ,V1 ,m1 , r/20); // jupiter
  for(t=tdibujo=0;t<20*T;t+=Deltat,tdibujo+=Deltat){
    
     
      // VELOCIDAD VERLET OPTIMIZADO ES DE CUARTO ORDEN, PERO EL PRECIO QUE HAY QUE PAGAR ES QUE SE CALCULAN DOS VECES LA FUERZA
    for(i=0;i<N;i++) Planeta[i].Mueva_r1(Deltat);
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<N;i++){Planeta[i].Mueva_V(Deltat);  Planeta[i].Mueva_r2(Deltat);}
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<N;i++){Planeta[i].Mueva_V(Deltat);  Planeta[i].Mueva_r1(Deltat);}
    
    cout << Planeta[1].Getx() << " " <<Planeta[1].Gety() << endl;
  }

  return 0;
}

