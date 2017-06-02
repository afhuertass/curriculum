#include <iostream>
#include <cmath>
#include "Vector.h"
using namespace std;



const double Deltat=0.001;  
const double chi=0.1931833250;
const double Um2chi = 1 -2*chi;
const double eps= 1.0;
const double r0 = 10;
const double m = 1;


class Cuerpo;
class Colisionador; // para el punto siguiente 

class Cuerpo{
private:
  vector3D r,V,F;  double m,R; //Posicion, Vel, Frueza, masa y Radio

public:
  void Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0);
  void BorreFuerza(void);
  void IncrementeFuerza();
  void Mueva_r1(double dt);
  void Mueva_V(double dt);
  void Mueva_r2(double dt);
  void Dibujese(void);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0){ //No es una funcion suelta sino que pertence a Cuerpo
  r.cargue(x0,y0,0);  V.cargue(Vx0,Vy0,0);
  m=m0, R=R0;
}
void Cuerpo::BorreFuerza(void){
  F.cargue(0,0,0);
}
void Cuerpo::IncrementeFuerza(void){
  // F+=F0; en esta funcion calculamos la fuerza debida al potencial.
  BorreFuerza();
  vector3D Fc;
  double nor = norma(r); // norma del vector de posicion
  double magnitudF = ((12*eps)/nor)*( pow(r0/nor ,12) - pow(r0/nor,6) );
  //cout << nor << endl;
  F = r*(magnitudF/nor); // la fuerza tiene direccion de r 
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


int main(void){

  Cuerpo molecula;
  double t, tdibujo;
  double T = 100; // areglos
  double ct;
  //void Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0);
  double kbT = 0.025;
  for ( kbT = 0.025-0.002 ; kbT < 0.025+0.002 ; kbT += 0.0001) { 
  molecula.Inicie(10,0,sqrt(2*kbT/m),0,m ,2.5 );
  double c = 0;
  double rprom = 0;
  for(t=tdibujo=0 ; t < T ; t+=Deltat, tdibujo+=Deltat){
    //verlet optimizado
    
    molecula.IncrementeFuerza();
    
    molecula.Mueva_V(Deltat); molecula.Mueva_r2(Deltat);
    molecula.IncrementeFuerza();
    molecula.Mueva_V(Deltat); molecula.Mueva_r1(Deltat);
    
    //cout << molecula.Getx() << " " << t << endl;
   //cout << t << " " << molecula.Getx() << endl;
       
    /*if( molecula.Getx() < 10 && c == 0 ){
      c++;
    }
    if ( c== 1 && molecula.Getx() > 10){
      cout << kbT <<" " <<t << endl; // periodo 
      break;
      } */ 
    
    // valor central de la posicion. 
    rprom += molecula.Getx();
    c+=1;
  }
  cout << kbT << " " << rprom/c << endl;
  
  }
  return 0.0;
}
