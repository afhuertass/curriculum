// Mi primer programa de DINAMICA MOLECULAR :)

#include <iostream>
#include <cmath>
#include "Random.h"
using namespace std;


 
const double tMAX = 400; // picosegundos
const double Deltat = 0.01; // picossegundos
const double L = 100; // Amstrongs
const double T = 298; // Kelvin
const double MASA = 22.89; // en u.m.a
const double e = 1; // carga
const double D = 0.132; // constante de difusion
const double kb = 0.826; // constate de bolztmann en esas unidades. 
const double Gamma = kb*T/(MASA*D);  
const double sigma = sqrt(2*D*Deltat);
const int N = 50;
class Cuerpo{
private: 
  double x,Vx,Fx, Fxpunto ,m,R, q; 
public:
  void Inicie(double x0, double Vx0, double m0, double R0 , double Q0);
  void CalculaFuerza(double Ex, double dt);      
  void Muevase(double dt , Crandom & ran2);
  void Dibujese(void);
  double Getx(void){return x;}; //Funcion inline
   
}; 



void Cuerpo::Inicie(double x0, double Vx0, double m0, double R0, double Q0 ){
  x=x0; Vx=Vx0; m=m0; R=R0;
  q = Q0;
  Fxpunto = 0;
  Fx = 0;
}

void Cuerpo::CalculaFuerza(double Ex, double dt){
  double Fxnew = q*Ex;
  Fxpunto = (Fxnew - Fx)/dt;
  Fx=q*Ex; // fuerza de campo E 
}

void Cuerpo::Muevase(double dt , Crandom & ran2){
  /*x+=Vx*dt;           y+=Vy*dt;
  Vx+=(Fx*dt)/m;     Vy+=(Fy*dt)/m;
  */
  double xnew = x + ( Fx + 0.5*Fxpunto*dt )*dt/(m*Gamma) + ran2.gauss(0, 2*D*dt );
  Vx = (xnew-x)/dt; x = xnew;
}

void Cuerpo::Dibujese(void){
  cout<<", "<<x<<"+"<<R<<"*cos(t),"<<0<<"+ "<<R<<"*sin(t)"<<endl;
}




//------FUNCIONES GLOBALES----------------------------

void InicieAnimacion(void){
  // cout<<"set terminal gif"<<endl;
  //cout<<"set output 'UnBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:110]"<<endl;
  cout<<"set yrange[-10:10]"<<endl;
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
  Cuerpo Balon[N];
  int i = 0;
  Crandom ran2(10);
  double t,tdibujo;
  double E, Jprom; 
  E = 0;
 
  //InicieAnimacion();
  //-------   (x0,y0,Vx0,Vy0,m0   ,R0);
  for(i = 0; i< N ; i++) Balon[i].Inicie(0 ,0 , MASA ,4 ,e);
  // for(t=tdibujo=0;t<tMAX;t+=Deltat,tdibujo+=Deltat){
  // if(tdibujo>tMAX/200){
  for(t=0;t<tMAX;t+=Deltat){
    //InicieCuadro();
    //for(i= 0 ; i < N ; i++)Balon[i].Dibujese();
    //TermineCuadro();
    for(i = 0; i< N ; i++ )cout  << t <<" "<<Balon[i].Getx()<<endl;
    for(i = 0; i< N ; i++) Balon[i].CalculaFuerza( E , Deltat);
    for(i = 0; i< N ; i++ ) Balon[i].Muevase(Deltat , ran2);
  }
  //}


  return 0;
  
}
