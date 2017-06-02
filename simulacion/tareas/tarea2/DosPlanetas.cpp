#include<iostream>
#include<cmath>
#include "Vector.h"

using namespace std;

const double Deltat=0.001;
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
  void DibujeseRelativo( Cuerpo sol, Cuerpo Jupiter );
  double ImprimirRotado(Cuerpo sol, Cuerpo Jupiter);
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
void Cuerpo::DibujeseRelativo(Cuerpo planeta1 , Cuerpo planeta2){
  // esta funcion recibe el sol y a jupiter, calcula el angulo entre los dos
  // y rota el presente cuerpo en dicho angulo
  vector3D dif = planeta2.GetR() - planeta1.GetR();
  double dx = planeta2.Getx() - planeta1.Getx();
  double dy = planeta2.Gety() - planeta1.Gety();
  double tan = dy/dx;
  double angle;
  if( dx == 0 && dy == 0){
    cout << ", " <<0 << "+" <<R<<"*cos(t) ,"<<0<<"+"<<R<<"*sin(t)";  
    return ;
  }
  if (dx < 0 ){
    tan = tan*-1;
    angle = atan(tan);
    angle = M_PI - angle;
      
  }else {
    angle = atan(tan);
  }
  vector3D dif2 = r - planeta1.GetR();
  double x = dif2.x()*cos(angle) + dif2.y()*sin(angle);
  double y = -dif2.x()*sin(angle) + dif2.y()*cos(angle);
  //
  cout << ", " <<x << "+" <<R<<"*cos(t) ,"<<y<<"+"<<R<<"*sin(t)";        
}
double Cuerpo::ImprimirRotado(Cuerpo planeta1 , Cuerpo planeta2){
  // adicionalmente retorna el valor de x 
  vector3D dif = planeta2.GetR() - planeta1.GetR();
  double dx = planeta2.Getx() - planeta1.Getx();
  double dy = planeta2.Gety() - planeta1.Gety();
  double tan = dy/dx;
  double angle;
  if( dx == 0 && dy == 0){
    //cout << ", " <<0 << "+" <<R<<"*cos(t) ,"<<0<<"+"<<R<<"*sin(t)";  
    cout << "0 0" << endl;
    return 0.0;
  }
  if (dx < 0 ){
    tan = tan*-1;
    angle = atan(tan);
    angle = M_PI - angle;
      
  }else {
    angle = atan(tan);
  }
  vector3D dif2 = r - planeta1.GetR();
  double x = dif2.x()*cos(angle) + dif2.y()*sin(angle);
  double y = -dif2.x()*sin(angle) + dif2.y()*cos(angle);
  //
  //cout << ", " <<x << "+" <<R<<"*cos(t) ,"<<y<<"+"<<R<<"*sin(t)";
  //cout << x << " " << y << endl; 
  return x;
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
  cout<<"set terminal png"<<endl;
  cout<< "set output 'res.png' " << endl;
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

  M=m0+m1;   omega=sqrt(G*M*pow(r,-3));  T=2*M_PI/omega;  
  r1=m0*r/M;  r0=r1-r;   V1=omega*r1;   V0=omega*r0;

  // InicieAnimacion();
  //void Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0);

  Planeta[0].Inicie(r0, 0 ,0 ,V0 ,m0 , r/10); // sol
  Planeta[1].Inicie(r1, 0 ,0 ,V1 ,m1 , r/20); // jupiter
  Planeta[2].Inicie(r1*cos(ang), r1*sin(ang) ,-V1*cos(ang2),V1*sin(ang2),mt , r/30 ); // el troyano 
  int k = 0;
  cout << "Periodo:" << T << endl;
  for(k=t=tdibujo=0;t<40*T;t+=Deltat,tdibujo+=Deltat , k++){
    
    if(k == 1000150000  ){
     InicieCuadro();
    for(i=0;i<N;i++) Planeta[i].Dibujese();
    for(i =0 ; i<N; i++) Planeta[i].DibujeseRelativo(Planeta[0], Planeta[1]);
    /* Planeta[1].DibujeseRelativo(Planeta[0] , Planeta[1]); // dibuje jupiter relativo al sol 
    Planeta[0].DibujeseRelativo(Planeta[0] , Planeta[1]);
    Planeta[2].DibujeseRelativo(Planeta[0], Planeta[1]); // troyano. 
    */
    TermineCuadro();
    tdibujo=0;
    }
    
    // cout<<Planeta[1].Getx()<<" "<<Planeta[1].Gety()<<endl;  //IMPRIME EN LA TERMINAL
    // VELOCIDAD VERLET OPTIMIZADO ES DE CUARTO ORDEN, PERO EL PRECIO QUE HAY QUE PAGAR ES QUE SE CALCULAN DOS VECES LA FUERZA
    for(i=0;i<N;i++) Planeta[i].Mueva_r1(Deltat);
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<N;i++){Planeta[i].Mueva_V(Deltat);  Planeta[i].Mueva_r2(Deltat);}
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<N;i++){Planeta[i].Mueva_V(Deltat);  Planeta[i].Mueva_r1(Deltat);}

    // obtener las coordenadas rotadas del troyano. 
    if ( k % 5000 == 0 ){
      double xr = Planeta[2].ImprimirRotado( Planeta[0] , Planeta[1]);
      //cout << t << " " << xr << endl;
    }
  }
  //cout << k << endl;
  return 0;
}

//  - experimientos para solucionar la tarea 
vector3D rotar(Cuerpo  planeta1 , Cuerpo  planeta2){
  // entonces el planeta1 es jupiter, el planeta2 es el sol.
  vector3D dif = planeta1.GetR() - planeta2.GetR();
  return dif;
  

}
