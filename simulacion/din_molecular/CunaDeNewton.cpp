// Mi primer pendulo :)

#include <iostream>
#include <cmath>
using namespace std;


const double GM=1;
const double g=979.9;
const double Deltat=0.000001;  
const double chi=0.193183325037836;


class Cuerpo{
private: 
  double theta,omega,N; //theta angulo, omega velocidad, N torque.
  double R,L,I,m,x0; //x0 posicion donde esta colgado 
  double xnew,ynew,Vxnew,Vynew; // variables nuevas de arranque
public:
  void Inicie(double theta0, double omega0,double m0, double R0,double L0,double x00);
  void CalculeTorqueBase(void);      
  void Mueva_r1(double dt);
  void Mueva_V(double dt);
  void Mueva_r2(double dt);
  void Dibujese(void);
  double Getx(void){return x0+L*sin(theta);}; //Funcion inline
  double Gety(void){return -L*cos(theta);};  
}; 



void Cuerpo::Inicie(double theta0, double omega0,double m0, double R0,double L0, double x00){
  theta=theta0;  omega=omega0;  R=R0; m=m0; L=L0, x0=x00;  I=m0*L*L;
}

//---------FUNCIONES NUEVAS --------------------------

void Cuerpo::Mueva_r1(double dt){  // se modifica el muevase
  theta+=omega*chi*dt;
   
 }
void Cuerpo::Mueva_V(double dt){  // se modifica el muevase
  omega+=(N*dt)/(2*I); 
}

void Cuerpo::Mueva_r2(double dt){  // se modifica el muevase
  theta+=omega*(1-2*chi)*dt;

}

void Cuerpo::CalculeTorqueBase(void){
  N=-L*m*g*sin(theta);   
}

//-----------------------------------------------

void Cuerpo::Dibujese(void){
  cout<<", "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+ "<<R<<"*sin(t)";
}




//------FUNCIONES GLOBALES----------------------------

void InicieAnimacion(void){
 
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:10]"<<endl;
  cout<<"set yrange[-14:0]"<<endl;
  // cout<<"set xrange["<<-L0<<":"<<L0<<"]"<<endl; cuando L0 es variable global
  //  cout<<"set yrange["<<-(L0+R0)*1.1<":0]"<<endl;
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
  double R0=1,M0=20,L0=10;
  double T;

 
  T=2*M_PI*sqrt(L0/g);
 
  
 
  InicieAnimacion();
 
  //-------   (        theta0, omega0, m0   ,R0, L0, x00);
  Pendulo.Inicie(15*M_PI/180, 0 ,M0 ,R0,L0 , 0);

  for(t=tdibujo=0;t<1.5*T*30;t+=Deltat,tdibujo+=Deltat){
    if(tdibujo>T/100){
    InicieCuadro();
    Pendulo.Dibujese();
    TermineCuadro();
    tdibujo=0;    }
    
   // cout<<Pendulo.Getx()<<" "<<Pendulo.Gety()<<endl;

   // se observa que se calcula dos veces la fuerza, la ventaja es que es de cuarto orden

      Pendulo.Mueva_r1(Deltat);
      Pendulo.CalculeTorqueBase();
      Pendulo.Mueva_V(Deltat);    Pendulo.Mueva_r2(Deltat);
      Pendulo.CalculeTorqueBase();
      Pendulo.Mueva_V(Deltat);   Pendulo.Mueva_r1(Deltat);
  }


  return 0;
  
}
