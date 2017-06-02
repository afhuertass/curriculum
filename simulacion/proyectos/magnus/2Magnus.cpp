#include<iostream>
#include<cmath>
#include "Vector.h"

using namespace std;
// Vertel
const double Deltat=0.001;
const double chi=0.193183325037836;
const double Um2chi=1-2*chi;

const double TMAX = 10;

// contantes del modelo
const double rho = 1.29;
const double CM = 10;
const double g = -10;
// rozamiento
const double b = 0.02; // friccion en el centro de masa
const double b2 = 0.0; // friccion en la rotacion 
const double eps = 0.001;

const int N=1;
//-------------------------------------------------------------
class Cuerpo;
class Colisionador;
//-------------------------------------------------------------

class Cuerpo{
private:
  vector3D r,V,F, omega , Tau;  double m,R; //Posicion, Vel, Frueza, masa y Radio
  double Area, I;
public:
  void Inicie(vector3D r0 , vector3D v0 ,vector3D omega0,double m0, double R0);
  void BorreFuerza(void);
  void IncrementeFuerza(vector3D F0);
  void IncrementeTorque(vector3D tau0);
  void Mueva_r1(double dt);
  void Mueva_V(double dt);
  void Mueva_Omega(double dt);
  void Mueva_r2(double dt);
  void Dibujese(void);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  double Getz(void){ return r.z() ;}; 
  void ImprimaCosas(void);
  friend class Colisionador;
};   // <-------- TIENE QUE COLOCARSE ESE ;!!!!!!!!!!

void Cuerpo::Inicie(vector3D r0 , vector3D v0 ,vector3D omega0,double m0, double R0){ //No es una funcion suelta sino que pertence a Cuerpo
  r.cargue(r0.x() , r0.y(),  r0.z() );  V.cargue(v0.x(), v0.y(), v0.z() );
  omega.cargue( omega0.x() , omega0.y(), omega0.z()  );
  m=m0 ; R=R0;
  Area = M_PI*pow(R,2);
  I = (2.0/5)*m*pow(R,2);
}
void Cuerpo::BorreFuerza(void){
  F.cargue(0,0,0);
  Tau.cargue(0,0,0);
}
void Cuerpo::IncrementeFuerza(vector3D F0){
  F+=F0;
}
void Cuerpo::IncrementeTorque(vector3D tau0){
  Tau+=tau0;
}
void Cuerpo::Mueva_r1(double dt){
  r+=V*(chi*dt); 
}
void Cuerpo::Mueva_V(double dt){
  V+=F*(dt/(2*m));
}
void Cuerpo::Mueva_Omega(double dt){
  omega+= Tau*(dt/(2*I));
}
void Cuerpo::Mueva_r2(double dt){
  r+=V*(Um2chi*dt);
}
void Cuerpo::Dibujese(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
void Cuerpo::ImprimaCosas(void){
  cout << omega.x() <<" "<<omega.y()<<" " <<omega.z() << endl;
}
//------------------------------------------------------------------
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Planeta);
  void AgregueFuerza(Cuerpo & Planeta1);
};

void Colisionador::CalculeFuerzas(Cuerpo * Planeta){
  int i,j;
  for(i=0;i<N;i++)
    Planeta[i].BorreFuerza();
  for(i=0;i<N;i++)
     AgregueFuerza(Planeta[i]);
}
void Colisionador::AgregueFuerza(Cuerpo & Planeta1){
  /*vector3D dr=Planeta1.r-Planeta2.r;
  double aux=G*Planeta1.m*Planeta2.m*pow(norma2(dr),-1.5);
  vector3D F2=dr*aux;
  Planeta2.IncrementeFuerza(F2);  Planeta1.IncrementeFuerza(F2*(-1));
  */
  vector3D aux = Planeta1.omega^Planeta1.V;
  vector3D gravedad;
  double magnitud = 0;
  gravedad.cargue(0,0 , g);
    aux = aux*(1/norma(aux));
  
  magnitud = 0.5*rho*Planeta1.Area*CM*norma2(Planeta1.V);
  aux = aux*(magnitud);
  //cout << norma(aux) << endl;

  Planeta1.IncrementeFuerza(aux);
  Planeta1.IncrementeFuerza(gravedad);
  
  // Fuerza de friccion  ;
  if( norma(Planeta1.V) > eps) {  // esto con el fin de evitar la fuerza de friccion cuando este quieto quizas sobre... 
  vector3D aux2  = (Planeta1.V )*(-1/norma(Planeta1.V));
  aux2 = aux2*b*norma2(Planeta1.V);
  Planeta1.IncrementeFuerza(aux2);

  }
  // una fuerza similar se propone para ir disminiyendo la velocidad angular
  if (norma(Planeta1.omega) > eps ){ 
  vector3D tau;
  tau.cargue(0,0,-1); // unitario en direccion menos z
  tau = tau*(Planeta1.R*Planeta1.omega.z()*b2); 
  Planeta1.IncrementeTorque(tau);
  }
  
}

//--------------------FUNCIONES GLOBALES-----------------------
void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl;
  //cout<<"set output 'MiPlaneta.gif'"<<endl;
  cout<<""<<endl;
  cout<<"unset key"<<endl;              // QUITAR LABEL
  cout<<"set xrange [-12:12]"<<endl;     // RANGO EN X
  cout<<"set yrange [-12:12]"<<endl;  
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
  double t,tdibujo , m0 ,R0 ;
  vector3D r0,v0, w0 ;
  int i = 0;
  m0 = 1; R0 = 0.01;
  r0.cargue(0,0,0);
  
  v0.cargue(50, 0, 10 );
  w0.cargue(0 ,0, 10 ) ;// Omega 
  Planeta[0].Inicie(r0, v0 , w0 , m0, R0);
  cout << "# Datos simulacion fuerza de magnus - Parametros" << endl;
  cout << "#velocidad inicial: " << norma(v0)<< endl;
  cout << "#  Omega_z " << w0.z() << ",  b:" << b << ", b2:" << b2 <<  endl;
  cout << "# Cm:" << CM << endl;
  for(t=tdibujo=0;t<TMAX;t+=Deltat,tdibujo+=Deltat){
    
    if(tdibujo>TMAX/500){
      /*InicieCuadro();
    for(i=0;i<N;i++) Planeta[i].Dibujese();
    TermineCuadro();
     tdibujo=0;
      */
    }
    
    // cout<<Planeta[1].Getx()<<" "<<Planeta[1].Gety()<<endl;  //IMPRIME EN LA TERMINAL
    // VELOCIDAD VERLET OPTIMIZADO ES DE CUARTO ORDEN, PERO EL PRECIO QUE HAY QUE PAGAR ES QUE SE CALCULAN DOS VECES LA FUERZA
    for(i=0;i<N;i++) Planeta[i].Mueva_r1(Deltat);
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<N;i++){Planeta[i].Mueva_V(Deltat); Planeta[i].Mueva_Omega(Deltat)  ; Planeta[i].Mueva_r2(Deltat);}
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<N;i++){Planeta[i].Mueva_V(Deltat); Planeta[i].Mueva_Omega(Deltat)  ;Planeta[i].Mueva_r1(Deltat);}

    cout << Planeta[0].Getx() << " " << Planeta[0].Gety() << " "<< Planeta[0].Getz() << endl ; 
    
    if( Planeta[0].Getz() <= 0) {  cout << t << endl ;break ; }
    // Planeta[0].ImprimaCosas();
  }

  return 0;
  
}
