//Mi Primer Programa de Dinamica Molecular
#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random.h"
using namespace std;

const double Deltat=0.01;
const double chi=0.193183325037836;
const double Um2chi=1-2*chi;

const double VEL0=10, OMEGA0=1.0;
const double K=50, Gamma=1.5, GammaTau=Gamma*0.6, Kcundall=K*0.6, mu=0.4;
const double g=1.0;

const double Lx=100,Ly=100;
const int Nx=5,Ny=5,N=Nx*Ny;

//---------------------------------------
class Cuerpo;
class Colisionador;
//---------------------------------------
class Cuerpo{
private:
  vector3D r,V,F; double m,R;  double theta,omega,tau,I;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double theta0,double omega0,double m0,double R0);
  void Mueva_r1(double dt);
  void Mueva_V(double dt);
  void Mueva_r2(double dt);
  void Dibujese(void);
  double Getx(void){return r.x();};//Funcion inline
  double Gety(void){return r.y();};//Funcion inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double theta0,double omega0,double m0,double R0){
  r.cargue(x0,y0,0); V.cargue(Vx0,Vy0,0);  theta=theta0; omega=omega0; I=2.0/5*m0*R0*R0;
  m=m0; R=R0;
}
void Cuerpo::Mueva_r1(double dt){
  r+=V*(chi*dt);  theta+=omega*chi*dt;
}
void Cuerpo::Mueva_V(double dt){
  V+=F*(dt/(2*m));  omega+=tau/I*dt/2;
}
void Cuerpo::Mueva_r2(double dt){
  r+=V*(Um2chi*dt);  theta+=omega*Um2chi*dt;
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
      <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t"; 
}
//---------------------------------------------
class Colisionador{
private:
  double dcontacto[N+4][N+4]; bool YaHabiaContacto[N+4][N+4];
public:
  void Inicie(void);
  void CalculeFuerzas(Cuerpo * Grano,double dt);
  void AgregueFuerza(Cuerpo & Grano1, Cuerpo & Grano2,int i,int j,double dt);
};
void Colisionador::Inicie(void){
  for(int i=0;i<N+4;i++)  
    for(int j=0;j<N+4;j++)
      {dcontacto[i][j]=0; YaHabiaContacto[i][j]=false;}
}
void Colisionador::CalculeFuerzas(Cuerpo * Grano,double dt){
  int i,j;
  for(i=0;i<N+4;i++)
    {Grano[i].F.cargue(0,0,0); Grano[i].tau=0;} 
  for(i=0;i<N;i++)
    {Grano[i].F.cargue(0,-Grano[i].m*g,0); Grano[i].tau=-GammaTau*Grano[i].omega;}
  for(i=0;i<N+4;i++)
    for(j=0;j<i&&j<N;j++)
      AgregueFuerza(Grano[i],Grano[j],i,j,dt);
}
void Colisionador::AgregueFuerza(Cuerpo & Grano1, Cuerpo & Grano2,int i,int j,double dt){
  vector3D r12,V12,Rw,Vc,n,t; double h,normar12,Fn,Ft,m1,m2,m12,R1,R2,omega1,omega2,Vn,Vt;
  r12=Grano1.r-Grano2.r; normar12=norma(r12);
  h=Grano1.R+Grano2.R-normar12; 
  if(h>0  && YaHabiaContacto[i][j]==false) YaHabiaContacto[i][j]=true;
  if(h<=0 && YaHabiaContacto[i][j]==true)  dcontacto[i][j]=0;

  if(h>0){
    //Geometría y dinámica del contacto
    m1=Grano1.m;   m2=Grano2.m;   m12=(m1*m2)/(m1+m2);
    R1=Grano1.m;   R2=Grano2.m;
    omega1=Grano1.omega;   omega2=Grano2.omega;

    V12=Grano1.V-Grano2.V;
    n=r12/normar12;  Rw.cargue(0,0,R1*omega1+R2*omega2);
    Vc=V12-(Rw^n); Vn=Vc*n; t.cargue(n.y(),-n.x(),0); Vt=Vc*t;

    //Fuerza de Hertz
    Fn=K*pow(h,1.5);
    //Disipacion
    Fn-=m12*sqrt(h)*Gamma*Vn; if(Fn<0) Fn=0;
    //Fuerza de fricción
    dcontacto[i][j]+=Vt*dt; Ft=Kcundall*dcontacto[i][j]; if (fabs(Ft)>mu*Fn) Ft=mu*Fn;
    
    //Cargar Resultados
    Grano1.F+=n*Fn;  Grano2.F-=n*Fn;
    Grano1.F-=t*Ft;  Grano2.F+=t*Ft;
    Grano1.tau-=Grano1.R*Ft;
    Grano2.tau-=Grano2.R*Ft;

  }
}
//----------- Funciones Globales --------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'MiBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-10:110]"<<endl;
  cout<<"set yrange [-10:110]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(double ytapa){
  cout<<"plot 100/7.0*t,0 , 0,100/7.0*t , 100,100/7.0*t , 100/7.0*t,"<<ytapa<<" ";
}
void TermineCuadro(void){
    cout<<endl;
}


int main(){
  Cuerpo Grano[N+4];
  Colisionador Newton;
  Crandom ran2(0);
  double t,tdibujo;
  int i,j,n;
  double dx=Lx/(Nx+1),dy=Ly/(Ny+1),theta,omega;
  double m0=0.1, R0=dx/4, T=40;

  //Las Paredes
  //-------------(x0      ,y0      ,Vx0,Vy0,theta0,omega0,m0,R0);
  Grano[N].Inicie( 10000+Lx,    Ly/2,0,0,0,0  ,1000*m0,10000); //pared east
  Grano[N+1].Inicie(   Lx/2,10000+Ly,0,0,0,0  ,1000*m0,10000); //pared north
  Grano[N+2].Inicie(-10000,     Ly/2,0,0,0,0  ,1000*m0,10000); //pared west
  Grano[N+3].Inicie(  Lx/2,   -10000,0,0,0,0  ,1000*m0,10000); //pared south
  //Los Granos
  for(n=0;n<N;n++){
    i=n%Nx; j=n/Nx; theta=2*M_PI*ran2.r(),omega=OMEGA0*(2*ran2.r()-1);
    //-------------(x0      ,y0      ,Vx0            ,Vy0            ,theta0,omega0,m0,R0);
    Grano[n].Inicie((i+1)*dx,(j+1)*dy,VEL0*cos(theta),VEL0*sin(theta),0     ,omega ,m0,R0);
  }
  
  InicieAnimacion();
  for(t=tdibujo=0;t<T;t+=Deltat,tdibujo+=Deltat){
    
    //    if(tdibujo>T/500){
    InicieCuadro(Grano[N+1].Gety()-10000);
    for(i=0;i<N;i++) Grano[i].Dibujese();
    TermineCuadro();
    tdibujo=0;
      //    }
    
    //cout<<t<<" "<<Grano[0].Gety()<<endl;

    for(i=0;i<N;i++) Grano[i].Mueva_r1(Deltat);
    Newton.CalculeFuerzas(Grano,0.5*Deltat);
    for(i=0;i<N;i++) {Grano[i].Mueva_V(Deltat); Grano[i].Mueva_r2(Deltat);}
    Newton.CalculeFuerzas(Grano,0.5*Deltat);
    for(i=0;i<N;i++) {Grano[i].Mueva_V(Deltat); Grano[i].Mueva_r1(Deltat);}
    
  }

  return 0;
}
