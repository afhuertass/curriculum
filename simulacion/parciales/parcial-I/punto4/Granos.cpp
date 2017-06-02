#include<iostream>
#include<cmath>
#include "Vector.h"
#include "Random.h"
using namespace std;

const double Deltat=0.01;
const double chi=0.193183325037836;
const double Um2chi=1-2*chi;
const double kbT = 0.5;
const double eps = 1.0;
const double r0 = 10;
const double m = 1; 

const double VEL0 = sqrt(2*kbT );
const int K = 100;

const double Lx = 100 , Ly = 100;
const int Nx = 5 , Ny = 5 , N = Nx*Ny;


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
  friend class Colisionador;
};   // <-------- TIENE QUE COLOCARSE ESE ; !!!!!!!!!!

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
  for(i=0;i<N+4;i++)
    Planeta[i].BorreFuerza();
  for(i=0;i<N+4;i++) 
    for(j=0;j<i && j< N;j++)
      AgregueFuerza(Planeta[i],Planeta[j]);
}
void Colisionador::AgregueFuerza(Cuerpo & Planeta1, Cuerpo & Planeta2){
  vector3D dr=Planeta1.r-Planeta2.r, Fn;
  double h, nor;
  nor = norma(dr);
  h = Planeta1.R + Planeta2.R - nor;
  /*if ( h > 0 ) {
    Fn = dr*(K*pow(h,1.5)/normadr);
    Planeta1.IncrementeFuerza(Fn) ; Planeta2.IncrementeFuerza(Fn*(-1));

    }*/
  // 
  double magnitud =  ((12*eps)/nor)*( pow(r0/nor ,12) - pow(r0/nor,6) );
  Fn = dr*(magnitud)/nor;
  Planeta1.IncrementeFuerza(Fn); Planeta2.IncrementeFuerza(Fn*(-1));
}

//--------------------FUNCIONES GLOBALES-----------------------
void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl;
  //cout<<"set output 'MiPlaneta.gif'"<<endl;
  cout<<""<<endl;
  cout<<"unset key"<<endl;              // QUITAR LABEL
  cout<<"set xrange [-10:110]"<<endl;     // RANGO EN X
  cout<<"set yrange [-10:110]"<<endl;  
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
  Crandom ran2(0);
  Cuerpo Planeta[N+4];
  Colisionador Newton;
  double t,tdibujo;
  double m0 = 0.1 , R0 ,T = 10 ;

  int i , n ,j;
  double dx = Lx/(Nx + 1)  ,dy = Ly/(Ny + 1) , theta;
  R0 = dx/4;
  // las paredes 
  /*Planeta[N].Inicie(10000 + Lx , Ly/2 , 0 ,0, 1000*m0, 10000 );
  Planeta[N+1].Inicie(Lx/2  , Ly + 10000 , 0 ,0, 1000*m0, 10000 );
  Planeta[N+2].Inicie(-10000 , Ly/2 , 0 ,0, 1000*m0 , 10000 );
  Planeta[N+3].Inicie( Lx/2 , -10000 , 0 ,0, 1000*m0, 10000 );
  */
  for (n = 0 ; n < N ; n++){
    
    i = n%Nx; j = n/Ny;
    theta = ran2.r()*2*M_PI;
    Planeta[n].Inicie( (i+1)*dx  , (j+1)*dy , VEL0*cos(theta) ,  VEL0*sin(theta), m0 , R0    );
    
  }

  
  InicieAnimacion();
 

  for(t=tdibujo=0;t<T;t+=Deltat,tdibujo+=Deltat){
    
  //if(tdibujo>T/500){
    InicieCuadro();
    for(i=0;i<N;i++) Planeta[i].Dibujese();
    TermineCuadro();
    //tdibujo=0;
    
    
    // cout<<Planeta[1].Getx()<<" "<<Planeta[1].Gety()<<endl;  //IMPRIME EN LA TERMINAL
    // VELOCIDAD VERLET OPTIMIZADO ES DE CUARTO ORDEN, PERO EL PRECIO QUE HAY QUE PAGAR ES QUE SE CALCULAN DOS VECES LA FUERZA
    for(i=0;i<N;i++) Planeta[i].Mueva_r1(Deltat);
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<N;i++){Planeta[i].Mueva_V(Deltat);  Planeta[i].Mueva_r2(Deltat);}
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<N;i++){Planeta[i].Mueva_V(Deltat);  Planeta[i].Mueva_r1(Deltat);}
  
  }
  return 0;
}
