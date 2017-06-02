// Difusion 1D 
#include <iostream>
#include <fstream>
#include <cmath>
// #include "Random64.h"
using namespace std;
// Automata difusion continua
const int Lx=200 ,Ly = 200 ; // tamaño cuadricula 

const double C =  0.5 ;// C < sqrt( 3/5)
const double C2 =C*C;
const double Um53C2 = 1- (5.0/3)*C2;
const double tau = 0.5;

const double A = 10;
const double lambda = 10;
const double omega =2*(2*M_PI*C)/lambda; 

class LatticeBolztmann{
private:
  double f[Lx][Ly][9], fnew[Lx][Ly][9]; // 
  int V[2][9]; 
  double w[9] ; // pesos de la funcion de equilibrio. 
public:
  LatticeBolztmann(void);
  void Inicie(void);
  void Muestre(char * nombreArchivo, int t);
  double rho(int ix, int iy, int t);
  double Jx(int ix, int iy);
  double Jy(int ix, int iy);
  double feq(double rho0, double Jx0 , double Jy0 , int i ); 
  void Colisione(int t);
  void Adveccione(void);
};

void  LatticeBolztmann::Adveccione(void){
  int i, ix, iy;
  for( ix = 0; ix <Lx ; ix++)
    for( iy = 0; iy <Ly ; iy++)
      for(i = 0 ; i < 9 ; i++)
	    f[ (ix + V[0][i]+Lx)%Lx ][ (iy + V[1][i] + Ly )%Ly ][i] = fnew[ix][iy][i];
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
      for(i =0 ; i < 9 ; i++)
	// evolucion ecuacion de bolztmann BGK
	fnew[ix][iy][i] = f[ix][iy][i] - 1.0/tau*( f[ix][iy][i] - feq(rho0,Jx0,Jy0,i) );  
    }
}

double LatticeBolztmann::feq(double rho0, double Jx0 , double Jy0 , int i ){
  
  if(i == 0)
    return rho0*Um53C2;
  else 
    return 3*w[i]*(C2*rho0 + V[0][i]*Jx0 + V[1][i]*Jy0 );
} 

double LatticeBolztmann::rho(int ix, int iy, int t){
  if( ix == Lx/2 && iy == Ly/2) 
    return  A*sin(omega*t);
  else { 
    int i; double suma;
    for(suma=0,i = 0; i < 9; i++)
      suma += f[ix][iy][i];
    return suma;
  }
}

double LatticeBolztmann::Jx(int ix, int iy){
  int i; double suma;
  for(suma=0,i = 0; i < 9; i++){
    suma += V[0][i]*f[ix][iy][i];
  }
  return suma;
}
double LatticeBolztmann::Jy(int ix, int iy){
  
  int i; double suma;
  for(suma=0,i = 0; i < 9; i++){
    suma += V[1][i]*f[ix][iy][i];
  }
  return suma;
}

LatticeBolztmann::LatticeBolztmann(void){
  // definir el valor de los pesos:
  w[0] = 4.0/9; w[1] = w[2] = w[3] = w[4] = 1.0/9;
  w[5] = w[6] = w[7] = w[8] = 1.0/36;
  
  // vectors de velocidad V[x,y][1...9]
  V[0][0] = 0;  V[0][1] = 1;   V[0][2] = 0; V[0][3] = -1;
  V[1][0] = 0;  V[1][1] = 0;   V[1][2] = 1 ; V[1][3] =0;

  V[0][4] = 0;  V[0][5] = 1;   V[0][6] = -1; V[0][7] = -1;
  V[1][4] = -1;  V[1][5] = 1;   V[1][6] = 1 ; V[1][7] = -1;

  V[0][8] = 1;
  V[1][8] = -1;
}
void LatticeBolztmann::Inicie(void){
  int i,ix,iy ;
  for( ix = 0 ; ix < Lx ; ix++)
    for( iy =0 ; iy< Ly ; iy++)
      for( i = 0 ; i < 9 ; i++)
	fnew[ix][iy][i] = f[ix][iy][i] = 0;
 
}

void LatticeBolztmann::Muestre(char *nombreArchivo, int t){
  ofstream MiArchivo(nombreArchivo);
  int ix,iy;
  for(ix = 50; ix < 150 ; ix++) {
    for( iy = 50 ; iy < 150 ; iy++)
      //cout << ix << " " << iy << endl;
      MiArchivo << ix << " "<< iy << " " << rho(ix,iy, t) <<endl;

  MiArchivo << endl;
  }
}

const int tmax=100;

int main(){
  LatticeBolztmann Ondas;
  int t;
  
  Ondas.Inicie();
  //Ondas.Muestre("Ondas.dat");
  for ( t = 0 ; t < 100 ; t++){
    
    Ondas.Adveccione();
    Ondas.Colisione(t);

  }
  Ondas.Muestre("Ondas2.dat", t);

  return 0;
}
