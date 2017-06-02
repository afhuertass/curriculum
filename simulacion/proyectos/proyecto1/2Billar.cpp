#include<iostream>
#include<cmath>
#include "Vector.h"

using namespace std;

const double Deltat=0.01;
const double chi=0.193183325037836;
const double Um2chi=1-2*chi;

const double VEL0 = 10;
const int K = 100;

const int N = 2;

const double mu = 0.00;
const double g = 9.8;
const double omega = 1;
const double nu = 0.3;

const double eps = 0.0001;
//-------------------------------------------------------------
class Cuerpo;
class Colisionador;
//-------------------------------------------------------------

class Cuerpo{
private:
  vector3D r,V,F;  double m,R; //Posicion, Vel, Frueza, masa y Radio
  
  bool colision;
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
  vector3D GetVelocidad(void) { return V;};
  bool GetColision(void){ return colision; };
  friend class Colisionador;
};   // <-------- TIENE QUE COLOCARSE ESE ; !!!!!!!!!!

void Cuerpo::Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0){ //No es una funcion suelta sino que pertence a Cuerpo
  r.cargue(x0,y0,0);  V.cargue(Vx0,Vy0,0);
  BorreFuerza();
  m=m0;  R=R0;
  colision = false;
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
  void FuerzaRozamiento(Cuerpo * Planeta);
};
void Colisionador::FuerzaRozamiento(Cuerpo * Planeta){
  double normaV = 0, magnitud = 0;
  vector3D Fr;
  for(int i = 0 ; i < 2 ; i++){
    //
    normaV = norma(Planeta[i].V);
    if (normaV < eps ) continue;
    if( Planeta[i].GetColision() ){
      // aplicamos fuerza de rozamiento, la fuerza de rozamiento
      magnitud = Planeta[i].m*g*mu;
      Fr = (Planeta[i].V)*(magnitud/normaV)*(-1);
      Planeta[i].IncrementeFuerza(Fr);
    }
  }

}
void Colisionador::CalculeFuerzas(Cuerpo * Planeta){
  int i,j;
  for(i=0;i<2;i++)
    Planeta[i].BorreFuerza();
  for(i=0;i<2;i++) 
    for(j=0;j<i && j< N;j++)
      AgregueFuerza(Planeta[i],Planeta[j]);
  FuerzaRozamiento(Planeta);
}
void Colisionador::AgregueFuerza(Cuerpo & Planeta1, Cuerpo & Planeta2){
  vector3D dr=Planeta1.r-Planeta2.r, Fn , Fr ;
  vector3D Fp ;
  double h, normadr , normaV;
  double normaFp , k2;
  normadr = norma(dr);
  h = Planeta1.R + Planeta2.R - normadr;
  if ( h > 0 ) {
    Fn = dr*(K*pow(h,1.5)/normadr);
    Planeta1.IncrementeFuerza(Fn) ; Planeta2.IncrementeFuerza(Fn*(-1));
    // vamos a agregar una fuerza perpendicular, constante mi
    k2 = -1*dr.x()/dr.y();
    Fp.cargue( 1 , k2 , 0);
    Fp = Fp*(1/norma(Fp)*nu*omega);
    //cout  << Fp*Fn << endl;
    Planeta2.IncrementeFuerza(Fp); 
    Planeta1.IncrementeFuerza(Fp*(-1));
    //Planeta1.IncrementeFuerza(Fn*(-10));
     if(Planeta1.colision == false ) {
      Planeta1.colision = true;
      Planeta2.colision = true;
    }

  }
  // fuerza de rozamiento
  
  //  Planeta2.IncrementeFuerza(F2);  Planeta1.IncrementeFuerza(F2*(-1));
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

double AnguloVectores( Cuerpo bola1 , Cuerpo bola2){
  double cosa = 0; // cosa del angulo
  double n1n2 = norma(bola1.GetVelocidad() )* norma(bola2.GetVelocidad() ) ;
  if(n1n2 == 0){
    return 0.0;
  }
  cosa = (bola1.GetVelocidad()*bola2.GetVelocidad())/n1n2;

  return acos(cosa);
  
}
int main(){
  
  Cuerpo Planeta[2];
  Colisionador Newton;
  double t,tdibujo;
  double m0 = 0.1 , R0 ,T = 50; 
  double v0 = 5;
  double beta = 0.99, y0 , b;
  double cosa = 0;
  double tc = 0;
  int i , n ,j , c = 0;

  R0 = 10;
  y0 = 50; 
  //cout << "#angulos -final simulacion -  mu= 0.00 , velocidad 5" << endl;
  for ( beta = 0.0 , c = 0; beta < 0.9 ; beta+= 0.0001) { 
  b = 2*R0*beta;
  Planeta[0].Inicie( 0, y0 , v0   , 0 , m0 , R0);
  Planeta[1].Inicie(30, y0-b, 0, 0, m0, R0);
  cosa = sqrt(  1 - (b*b)/(4*R0*R0) );
  //cout << b << endl;
  //InicieAnimacion();
  c = 0; // 
  for(t=tdibujo=0;t<T;t+=Deltat,tdibujo+=Deltat){
    /*InicieCuadro();
    Planeta[0].Dibujese();
    Planeta[1].Dibujese();
    TermineCuadro();
    */
    // VELOCIDAD VERLET OPTIMIZADO ES DE CUARTO ORDEN, PERO EL PRECIO QUE HAY QUE PAGAR ES QUE SE CALCULAN DOS VECES LA FUERZA
    for(i=0;i<2;i++) Planeta[i].Mueva_r1(Deltat);
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<2;i++){Planeta[i].Mueva_V(Deltat);  Planeta[i].Mueva_r2(Deltat);}
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<2;i++){Planeta[i].Mueva_V(Deltat);  Planeta[i].Mueva_r1(Deltat);}

    if( Planeta[0].GetColision() ){
     
      //break;
    } 
  }
  //
  cout << beta << " " << AnguloVectores(Planeta[0], Planeta[1])*(180.0/M_PI) << endl ;
    
   }
  
  return 0;
}