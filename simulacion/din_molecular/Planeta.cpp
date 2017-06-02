// Mi primer programa de DINAMICA MOLECULAR :)

#include <iostream>
#include <cmath>
using namespace std;


const double GM=1;
const double Deltat=0.1;  
//const double M_PI= ;

class Cuerpo{
private: 
  double x,y,Vy,Vx,Fx,Fy,m,R; 
public:
  void Inicie(double x0,double y0, double Vx0, double Vy0,double m0, double R0);
  void CalculaFuerza(void);      
  void Muevase(double dt);
  void Dibujese(void);
  double Getx(void){return x;}; //Funcion inline
  double Gety(void){return y;};  
}; 



void Cuerpo::Inicie(double x0,double y0, double Vx0, double Vy0,double m0, double R0){
  x=x0; y=y0;  Vy=Vy0; Vx=Vx0; m=m0; R=R0;
}

void Cuerpo::CalculaFuerza(void){
  double aux=-GM*m*pow(x*x+y*y,-1.5);
  Fx=aux*x;
  Fy=aux*y;
}

void Cuerpo::Muevase(double dt){
  x+=Vx*dt;           y+=Vy*dt;
  Vx+=(Fx*dt)/m;     Vy+=(Fy*dt)/m;
}

void Cuerpo::Dibujese(void){
  cout<<", "<<x<<"+"<<R<<"*cos(t),"<<y<<"+ "<<R<<"*sin(t)"<<endl;
}




//------FUNCIONES GLOBALES----------------------------

void InicieAnimacion(void){
 
  cout<<"unset key"<<endl;
  cout<<"set xrange[-12:12]"<<endl;
  cout<<"set yrange[-12:12]"<<endl;
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
  Cuerpo Balon;
  double t,tdibujo;
  double r=10,m=1;
  double omega,T,V;

  omega=sqrt(GM*pow(r,-3));
  T=2*M_PI/omega;
  V=omega*r;
  
 
  //InicieAnimacion();
  //-------   (x0,y0,Vx0,Vy0,m0   ,R0);
  Balon.Inicie(r ,0 ,0 ,V ,m , r/20);

  for(t=tdibujo=0;t<3*T;t+=Deltat,tdibujo+=Deltat){
    /*
    if(tdibujo>T/100){
    InicieCuadro();
    Balon.Dibujese();
    TermineCuadro();
  //  
  */
  cout<<Balon.Getx()<<" "<<Balon.Gety()<<endl;

   tdibujo=0;    

 Balon.CalculaFuerza();
       Balon.Muevase(Deltat);
  }


  return 0;
  
}
