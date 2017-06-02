// Difusion 1D 
#include <iostream>
#include <fstream>
#include <cmath>
// #include "Random64.h"
using namespace std;
// Automata difusion continua
const int Lx=100 ,Ly = 100 ; // tamaño cuadricula 

const double C =  0.5 ;// C < sqrt( 3/5)
const double C2 =C*C;
const double Um3Umw0  = 1-2*C2;
const double tau = 0.5;

const double A = 10;
const double lambda = 10;
const double omega =(2*M_PI*C)/lambda; 

const int ixc = 30 , iyc = Ly/2 , R = 50;
enum tipoCelda { normal,  abso, fuente, espejo, lente};
class LatticeBolztmann{
private:
  double f[Lx][Ly][5], fnew[Lx][Ly][5]; // 
  int V[2][5]; 
  double w[5] ; // pesos de la funcion de equilibrio. 
  tipoCelda CeldaTipo[Lx][Ly];  // cero = normal , 1 = absorbente,  2 = fuente , 3 = espejo 

public:
  LatticeBolztmann(void);
  void Inicie(void);
  void Muestre(char * nombreArchivo, int t);
  void Muestre2(int t);
  double rho(int ix, int iy, int t);
  double Jx(int ix, int iy);
  double Jy(int ix, int iy);
  double feq(double rho0, double Jx0 , double Jy0 , int i , double vel ); 
  void Colisione(int t);
  void Adveccione(void);
  void ConstruyaGeometria(void);
  double Velocidad(int ix, int iy );
  };
void LatticeBolztmann::ConstruyaGeometria(void){
  int ix, iy;
  for( ix = 0 ; ix < Lx ; ix++)
    for(iy = 0 ; iy < Ly ; iy++)
      CeldaTipo[ix][iy] = normal;
  ix = 0 ; for(iy = 0 ; iy < Ly ; iy++) CeldaTipo[ix][iy] = abso;
  /*ix = 0 ; for(iy = 0 ; iy < Ly ; iy++) CeldaTipo[ix][iy] = abso;
  iy = 0 ; for(ix = 0 ; ix < Lx ; ix++) CeldaTipo[ix][iy] = abso;
  
  iy = Ly-1; for(ix = 0 ; ix < Lx ; ix++) CeldaTipo[ix][iy] = abso;
  
  ix = Lx-1 ; for(iy = 0 ; iy < Ly ; iy++) CeldaTipo[ix][iy] = abso;
  
  */
 // rendijas 
 /*ix =50;  for ( iy =0 ; iy < Ly/2-5 ; iy++ ) CeldaTipo[ix][iy] = abso;
 for(iy = Ly/2+5 ; iy < Ly ; iy++) CeldaTipo[ix][iy] = abso;
 int ancho = 50, iy0 = 0;
 //for(iy = iy0 ; iy < iy0 + ancho; iy++) CeldaTipo[ix][iy] == normal;
 iy0 = 50;
 //for(iy = iy0 ; iy < iy0 + ancho; iy++) CeldaTipo[ix][iy] == normal;

 */
 ix = 1 ;
 for( iy = 0 ; iy < Ly; iy++)
   CeldaTipo[ix][iy] = fuente;
 //CeldaTipo[1][Ly-3] = fuente;
 
  

  ix = 0; iy = 0;
  for( ix = 0; ix < Lx ; ix++)
    for( iy = 0 ; iy < Ly ; iy++)
	   if( ix > 30 && ix < Lx) // entre 100 y 200
	     if ( pow( ix - ixc,2)+pow(iy-iyc,2) > pow(R,2)) // que esten fuera del circulo
	       if ( ix > ixc) // y a la derecha
		 CeldaTipo[ix][iy] = espejo; // son espejos
		 
}
double LatticeBolztmann::Velocidad(int ix, int iy ){
  /*if ( ix < 50 ){ imprementacion sin suavizado
    return  0.5;
  }else {
    return 0.5/2; 
  }
  */
  /*  double n = 0.5*tanh(ix - 50 )+1.5; // suavizado
  return 0.5/n;
  */
  // ecuacion de la recta que forma 20 grados con la vertical 0.36*iy+50
  //double n = 0.5*tanh(ix - ( 1.73*iy+10 ) )+1.5; // suavizado
  //return 0.5/n; // 
  return 0.5;
}

void  LatticeBolztmann::Adveccione(void){
  int i, ix, iy;
  for( ix = 0; ix <Lx ; ix++)
    for( iy = 0; iy <Ly ; iy++)
      for(i = 0 ; i < 5 ; i++) { 
	if( CeldaTipo[ix][iy] == abso) continue; // las absorbentes no adveccionan?
	f[ (ix + V[0][i]+Lx)%Lx ][ (iy + V[1][i] + Ly )%Ly ][i] = fnew[ix][iy][i];
      }
}

void LatticeBolztmann::Colisione(int t){
  int ix, iy,i;
  double rho0, Jx0, Jy0, vel;
  for(ix =0 ; ix < Lx ; ix++ )  //para cada celda
    for(iy =0 ; iy < Ly ; iy++) {
      rho0 = rho(ix,iy,t);
      Jx0 = Jx(ix,iy); 
      Jy0 = Jy(ix,iy);
      vel = Velocidad(ix,iy);
      // para cada direccion
      if(CeldaTipo[ix][iy] == espejo) {
	// aqui implentar la regla de los espejos
	double aux = 0;
	aux = f[ix][iy][1];
	f[ix][iy][1] = f[ix][iy][3]; f[ix][iy][3] = aux;
	aux = f[ix][iy][2];
	f[ix][iy][2] = f[ix][iy][4]; f[ix][iy][4] = aux;
	
	continue ;
      }
      for(i =0 ; i < 5 ; i++)
	// evolucion ecuacion de bolztmann BGK
	fnew[ix][iy][i] = f[ix][iy][i] - 1.0/tau*( f[ix][iy][i] - feq(rho0,Jx0,Jy0,i , vel) );  
    }
}

double LatticeBolztmann::feq(double rho0, double Jx0 , double Jy0 , int i, double vel ){
  double vel2 = vel*vel;
  if(i == 0)
    return rho0*(1-2*vel2 ); //  rho0*Um3Umw0;
  else 
    return 3*w[i]*(vel2*rho0 + V[0][i]*Jx0 + V[1][i]*Jy0 );
} 

double LatticeBolztmann::rho(int ix, int iy, int t){
  if(CeldaTipo[ix][iy] == fuente) 
    return  A*sin(omega*t);
  if(CeldaTipo[ix][iy] == abso)
    return 0.0;
  if(CeldaTipo[ix][iy] == normal)  { 
    int i; double suma;
    for(suma=0,i = 0; i < 5; i++)
      suma += f[ix][iy][i];
    return suma;
  }
}

double LatticeBolztmann::Jx(int ix, int iy){
  int i; double suma;
  for(suma=0,i = 0; i < 5; i++){
    suma += V[0][i]*f[ix][iy][i];
  }
  return suma;
}
double LatticeBolztmann::Jy(int ix, int iy){
  
  int i; double suma;
  for(suma=0,i = 0; i < 5; i++){
    suma += V[1][i]*f[ix][iy][i];
  }
  return suma;
}

LatticeBolztmann::LatticeBolztmann(void){
  // definir el valor de los pesos:
  w[0] = 1.0/3; w[1] = w[2] = w[3] = w[4] = 1.0/6;
  
  
  // vectors de velocidad V[x,y][1...9]
  V[0][0] = 0;  V[0][1] = 1;   V[0][2] = 0; V[0][3] = -1;
  V[1][0] = 0;  V[1][1] = 0;   V[1][2] = 1 ; V[1][3] =0;

  V[0][4] = 0;
  V[1][4] = -1;
}
void LatticeBolztmann::Inicie(void){
  int i,ix,iy ;
  for( ix = 0 ; ix < Lx ; ix++)
    for( iy =0 ; iy< Ly ; iy++)
      for( i = 0 ; i < 5 ; i++)
	fnew[ix][iy][i] = f[ix][iy][i] = 0;
 
  ConstruyaGeometria();
}

void LatticeBolztmann::Muestre(char *nombreArchivo, int t){
  ofstream MiArchivo(nombreArchivo);
  int ix,iy;
  for(ix = 0; ix < 100 ; ix+= 1  ) {
    for( iy = 0 ; iy < 100 ; iy+= 1)
      //cout << ix << " " << iy << endl;
      MiArchivo << ix << " "<< iy << " " << rho(ix,iy, t) <<endl;

  MiArchivo << endl;
  }
}
void LatticeBolztmann::Muestre2(int t){
  // obtener los 
  int ix = 0,iy = 0  ;
  ofstream MiArchivo("tipos.dat");
   for(ix = 0; ix < Lx ; ix+= 1  ) {
     for( iy = 0 ; iy < Ly ; iy+= 1) { 
      int a = 0;
      if(CeldaTipo[ix][iy] == espejo)
	a = 1;
      MiArchivo << ix <<" " << iy << " "  <<  a << endl;
     }
     MiArchivo <<endl;
   }
}
const int tmax=1000;

int main(){
  LatticeBolztmann Ondas;
  int t;
  
  Ondas.Inicie();
  //Ondas.Muestre("Ondas.dat");
  for ( t = 0 ; t < tmax ; t++){
    
    Ondas.Adveccione();
    Ondas.Colisione(t);
   
  }
  Ondas.Muestre("2DQ5-espejo.dat", t);
  Ondas.Muestre2(t);
  return 0;
}
