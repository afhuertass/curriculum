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
  double feq(double rho0, double Jx0 , double Jy0 , int i ); 
  void Colisione(int t);
  void Adveccione(void);
  void ConstruyaGeometria(void);
  };
void LatticeBolztmann::ConstruyaGeometria(void){
  int ix, iy;
  for( ix = 0 ; ix < Lx ; ix++)
    for(iy = 0 ; iy < Ly ; iy++)
      CeldaTipo[ix][iy] = normal;
 
  iy = 0 ; for(ix = 0 ; ix < Lx ; ix++) CeldaTipo[ix][iy] = abso;
  iy = Lx-1; for(ix = 0 ; ix < Lx ; ix++) CeldaTipo[ix][iy] = abso;
  
 ix = 0 ; for(iy = 0 ; iy < Ly ; iy++) CeldaTipo[ix][iy] = abso;
 ix = Lx-1 ; for(iy = 0 ; iy < Ly ; iy++) CeldaTipo[ix][iy] = abso;
 
 // rendijas 
 

 //CeldaTipo[1][0] = fuente;
 CeldaTipo[1][Ly-1] = fuente;
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
  double rho0, Jx0, Jy0;
  for(ix =0 ; ix < Lx ; ix++ )  //para cada celda
    for(iy =0 ; iy < Ly ; iy++) {
      rho0 = rho(ix,iy,t);
      Jx0 = Jx(ix,iy); 
      Jy0 = Jy(ix,iy);
      // para cada direccion
      for(i =0 ; i < 5 ; i++)
	// evolucion ecuacion de bolztmann BGK
	fnew[ix][iy][i] = f[ix][iy][i] - 1.0/tau*( f[ix][iy][i] - feq(rho0,Jx0,Jy0,i) );  
    }
}

double LatticeBolztmann::feq(double rho0, double Jx0 , double Jy0 , int i ){
  
  if(i == 0)
    return rho0*Um3Umw0;
  else 
    return 3*w[i]*(C2*rho0 + V[0][i]*Jx0 + V[1][i]*Jy0 );
} 

double LatticeBolztmann::rho(int ix, int iy, int t){
  if(CeldaTipo[ix][iy] == fuente) 
    {
      if( iy > Ly/2 ) {
      return  A*sin(omega*t);
      } else { // una de las fuentes esta desfazada
	return A*sin(omega*t+M_PI/2);
      }
    }
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
  int ix = 50,iy  ;
  for( iy = 0 ; iy < Ly ; iy++){
    cout << iy << " " << pow(rho(ix , iy , t),2 ) << endl;
  }

}
const int tmax=100;

int main(){
  LatticeBolztmann Ondas;
  int t;
  
  Ondas.Inicie();
  //Ondas.Muestre("Ondas.dat");
  for ( t = 0 ; t < 200 ; t++){
    
    Ondas.Adveccione();
    Ondas.Colisione(t);
        
  }
  //Ondas.Muestre("2DQ5-desfase.dat", t);
  Ondas.Muestre2(t);
  return 0;
}
