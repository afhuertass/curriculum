// Mi primer programa de DINAMICA MOLECULAR :)

#include <iostream>
#include <cmath>
using namespace std;


const double g=9.8;
const  double Deltat=0.001;  

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
  Fx=0; Fy=-m*g;
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
  // cout<<"set terminal gif"<<endl;
  //cout<<"set output 'UnBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[0:300]"<<endl;
  cout<<"set yrange[-100:100]"<<endl;
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
   
 
  InicieAnimacion();
  //-------   (x0,y0,Vx0,Vy0,m0   ,R0);
  Balon.Inicie(0 ,0 ,30 ,40 ,0.453,5);
  // for(t=tdibujo=0;t<tMAX;t+=Deltat,tdibujo+=Deltat){
  // if(tdibujo>tMAX/200){
  for(t=0;t<10;t+=Deltat){
    InicieCuadro();
    Balon.Dibujese();
    TermineCuadro();
    cout<<Balon.Getx()<<" "<<Balon.Gety()<<endl;
    Balon.CalculaFuerza();
    Balon.Muevase(Deltat);
  }
  //}


  return 0;
  
}
