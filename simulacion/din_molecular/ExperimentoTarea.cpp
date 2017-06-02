// Mi primer programa de DINAMICA MOLECULAR :)

#include <iostream>
#include <cmath>
using namespace std;


const double GM=1;
const double Deltat=0.01;  
const double chi=0.1931833250;

class Cuerpo{
private: 
  double x,y,xold , yold, Vy,Vx,Fx,Fy,m,R;
  
public:
  Cuerpo ();
  Cuerpo (double,double,double,double,double,double );
  void Inicie(double x0,double y0, double Vx0, double Vy0,double m0, double R0);
  void CalculaFuerza(void);      
  void Mueva_r1(double dt);
  void Mueva_r2(double dt);
  void Mueva_v(double dt);
  void Dibujese(void);
  double Getx(void){return x;}; //Funcion inline
  double Gety(void){return y;};  
}; 
Cuerpo::Cuerpo(){

}
Cuerpo::Cuerpo(double x0,double y0 , double Vx0 , double Vy0, double m0, double R0){
  this->Inicie(x0,y0,Vx0,Vy0,m0, R0);
}

void Cuerpo::Inicie(double x0,double y0, double Vx0, double Vy0,double m0, double R0){
  x=x0; y=y0;  Vy=Vy0; Vx=Vx0; m=m0; R=R0;
}

void Cuerpo::CalculaFuerza(void){
  double aux=-GM*m*pow(x*x+y*y,-1.5);
  Fx=aux*x;
  Fy=aux*y;
}

void Cuerpo::Mueva_r1(double dt){
  x+= Vx*chi*dt;y+= Vy*chi*dt;
}
void Cuerpo::Mueva_v(double dt){
  Vx += Fx*dt/(2*m); Vy += Fy*dt/(2*m);
}
void Cuerpo::Mueva_r2(double dt){
  x+= Vx*(1-2*chi)*dt;
  y+= Vy*(1-2*chi)*dt;
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
  Cuerpo * planetas[2] ; // un arreglo de planetas 
  omega=sqrt(GM*pow(r,-3));
  T=2*M_PI/omega;
  V=omega*r;
  //  
  planetas[0] = new Cuerpo(r,0,0,V,m,r/20); // primer planeta
  planetas[1] = new Cuerpo(r,0,0,V,m*20,r); // segundo planeta 
  // ahora necesitamos.... un metododo, que reciba un Cuerpo 
  // y con la masa y coordenadas calcule la fuerza 
  /*InicieAnimacion();
  //-------   (x0,y0,Vx0,Vy0,m0   ,R0);
  */
  Balon.Inicie(r ,0 ,0 ,V ,m , r/20);
 
  for(t=tdibujo=0;t<3000*T;t+=Deltat,tdibujo+=Deltat){
    
    cout<<Balon.Getx()<<" "<<Balon.Gety()<<endl;
    
    Balon.Mueva_r1(Deltat);
    Balon.CalculaFuerza();
    Balon.Mueva_v(Deltat);   Balon.Mueva_r2(Deltat);
    Balon.CalculaFuerza();
    Balon.Mueva_v(Deltat); Balon.Mueva_r1(Deltat);

   
  }
  

  return 0;
  
}
