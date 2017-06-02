// Difusion 1D 
#include <iostream>
#include <fstream>
#include <cmath>
// #include "Random64.h"
using namespace std;
// Automata difusion continua
const int Lx=100 ,Ly = 100 ; // tamaño cuadricula 
const int Lz = 100;

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
  double f[Lx][Ly][Lz][15], fnew[Lx][Ly][Lz][15]; // 
  int V[3][15]; 
  double w[15] ; // pesos de la funcion de equilibrio. 
  tipoCelda CeldaTipo[Lx][Ly][Lz];  // cero = normal , 1 = absorbente,  2 = fuente , 3 = espejo 

public:
  LatticeBolztmann(void);
  void Inicie(void);
  void Muestre(char * nombreArchivo, int t);
  void Muestre2(int t);
  double rho(int ix, int iy,int iz, int t);
  double Jx(int ix, int iy, int iz);
  double Jy(int ix, int iy ,int iz);
  double Jz(int ix,int iy , int iz);
  double feq(double rho0, double Jx0 , double Jy0 , double Jz0 , int i ); 
  void Colisione(int t);
  void Adveccione(void);
  void ConstruyaGeometria(void);
  };
void LatticeBolztmann::ConstruyaGeometria(void){
  int ix, iy, iz;
  for( ix = 0 ; ix < Lx ; ix++)
    for(iy = 0 ; iy < Ly ; iy++)
      for ( iz = 0 ; iz < Lz ; iz++)
	CeldaTipo[ix][iy][iz] = normal;
  /*
  iy = 0 ; for(ix = 0 ; ix < Lx ; ix++) CeldaTipo[ix][iy] = abso;
  iy = Lx-1; for(ix = 0 ; ix < Lx ; ix++) CeldaTipo[ix][iy] = abso;
  
 ix = 0 ; for(iy = 0 ; iy < Ly ; iy++) CeldaTipo[ix][iy] = abso;
 ix = Lx-1 ; for(iy = 0 ; iy < Ly ; iy++) CeldaTipo[ix][iy] = abso;
 
 // rendijas 
 
 */
  ix = 0; 
  for( iy = 0 ; iy < Ly ; iy++)
    for(iz = 0 ; iz < Lz ; iz++)
      CeldaTipo[ix][iy][iz] = abso;
  
  ix = Lx-1; 
  for( iy = 0 ; iy < Ly ; iy++)
    for(iz = 0 ; iz < Lz ; iz++)
      CeldaTipo[ix][iy][iz] = abso;
  
  iy = 0; 
  for( ix = 0 ; ix < Lx ; ix++)
    for(iz = 0 ; iz < Lz ; iz++)
      CeldaTipo[ix][iy][iz] = abso;
  
  iy = Ly-1; 
  for( ix = 0 ; ix < Lx ; ix++)
    for(iz = 0 ; iz < Lz ; iz++)
      CeldaTipo[ix][iy][iz] = abso;
  
  iz = 0; 
  for( ix = 0 ; ix < Lx ; ix++)
    for(iy = 0 ; iy < Ly ; iy++)
      CeldaTipo[ix][iy][iz] = abso;

  iz = Lz-1; 
  for( ix = 0 ; ix < Lx ; ix++)
    for(iy = 0 ; iy < Ly ; iy++)
      CeldaTipo[ix][iy][iz] = abso;
  
  
  CeldaTipo[Lx/2][0][Lz/2] = fuente;
  //CeldaTipo[1][Ly-1][0] = fuente;

 }

void  LatticeBolztmann::Adveccione(void){
  int i, ix, iy, iz ;
  for( ix = 0; ix <Lx ; ix++)
    for( iy = 0; iy <Ly ; iy++)
      for( iz = 0 ; iz < Lz ; iz++)
	for(i = 0 ; i < 15 ; i++) { 
	  if( CeldaTipo[ix][iy][iz] == abso) continue; // las absorbentes no adveccionan?
	  f[ (ix + V[0][i]+Lx)%Lx ][ (iy + V[1][i] + Ly )%Ly ][ (iz + V[2][i]+Lz)%Lz ][i] = fnew[ix][iy][iz][i];
	}
}

void LatticeBolztmann::Colisione(int t){
  int ix, iy, iz,i;
  double rho0, Jx0, Jy0 , Jz0;
  for(ix =0 ; ix < Lx ; ix++ )  //para cada celda
    for(iy =0 ; iy < Ly ; iy++) 
      for( iz = 0 ; iz < Lz ; iz++)
	{
	  rho0 = rho(ix,iy,iz,t);
	  Jx0 = Jx(ix,iy,iz); 
	  Jy0 = Jy(ix,iy,iz);
	  Jz0 = Jz(ix,iy,iz);
	  // para cada direccion
	  for(i =0 ; i < 15 ; i++)
	// evolucion ecuacion de bolztmann BGK
	    fnew[ix][iy][iz][i] = f[ix][iy][iz][i] - 1.0/tau*( f[ix][iy][iz][i] - feq(rho0,Jx0,Jy0,Jz0,i) );  
	}
}

double LatticeBolztmann::feq(double rho0, double Jx0 , double Jy0, double Jz0 , int i ){
  
  if(i == 0)
    return 9.0/2*w[0]*rho0*(1-7.0/3*C2);  //rho0*Um3Umw0;
  else 
    return 3*w[i]*(C2*rho0 + V[0][i]*Jx0 + V[1][i]*Jy0 + V[2][i]*Jz0 );
} 

double LatticeBolztmann::rho(int ix, int iy,int iz, int t){
  if(CeldaTipo[ix][iy][iz] == fuente) 
    return  A*sin(omega*t);
  if(CeldaTipo[ix][iy][iz] == abso)
    return 0.0;
  if(CeldaTipo[ix][iy][iz] == normal)  { 
    int i; double suma;
    for(suma=0,i = 0; i < 15; i++)
      suma += f[ix][iy][iz][i];
    return suma;
  }
}

double LatticeBolztmann::Jx(int ix, int iy, int iz){
  int i; double suma;
  for(suma=0,i = 0; i < 15; i++){
    suma += V[0][i]*f[ix][iy][iz][i];
  }
  return suma;
}
double LatticeBolztmann::Jy(int ix, int iy, int iz){
  
  int i; double suma;
  for(suma=0,i = 0; i < 15; i++){
    suma += V[1][i]*f[ix][iy][iz][i];
  }
  return suma;
}
double LatticeBolztmann::Jz(int ix,int iy , int iz){
 int i; double suma;
  for(suma=0,i = 0; i < 15; i++){
    suma += V[2][i]*f[ix][iy][iz][i];
  }
  return suma;
  
}

LatticeBolztmann::LatticeBolztmann(void){
  // definir el valor de los pesos:
  // w[0] = 1.0/3; w[1] = w[2] = w[3] = w[4] = 1.0/6;
   w[0] = 2.0/9;
  w[1] = w[2] = w[3] = w[4] = w[5] = w[6] =1.0/9;
  w[7] = w[8] = w[9] = w[10] = w[11] = w[12] = w[13] = w[14] = 1.0/72;
 
  
  V[0][0] = 0 ; V[1][0] = 0 ; V[2][0] = 0; 
  V[0][1] = 1 ; V[1][1] = 0 ; V[2][1]= 0;
  V[0][2] = 0 ; V[1][2] = 1 ; V[2][2]= 0;
  V[0][3] = -1 ; V[1][3] = 0 ; V[2][3]= 0;
  V[0][4] = 0 ; V[1][4] = -1 ; V[2][4]= 0;
  V[0][5] = 0 ; V[1][5] = 0 ; V[2][5]= 1;
  V[0][6] = 0 ; V[1][6] = 0 ; V[2][6]= -1;

  V[0][7] = 1 ; V[1][7] = 1 ; V[2][7]= 1;
  V[0][8] = -1 ; V[1][8] = 1 ; V[2][8]= 1;
  V[0][9] = -1 ; V[1][9] = -1 ; V[2][9]= 1;
  V[0][10] = 1 ; V[1][10] = -1 ; V[2][10]= 1;
  V[0][11] = 1 ; V[1][11] = 1 ; V[2][11] =-1;
 
  V[0][12] = -1 ; V[1][12] = 1 ; V[2][12]= -1;
  V[0][13] = -1 ; V[1][13] = -1 ; V[2][13]= -1;
  V[0][14] = 1 ; V[1][14] = -1 ; V[2][14]= -1;
  
}
void LatticeBolztmann::Inicie(void){
  int i,ix,iy,iz ;
  for( ix = 0 ; ix < Lx ; ix++)
    for( iy =0 ; iy< Ly ; iy++)
      for(iz = 0 ; iz < Lz ; iz++)
	for( i = 0 ; i < 15 ; i++)
	  fnew[ix][iy][iz][i] = f[ix][iy][iz][i] = 0;
 
  ConstruyaGeometria();
}

void LatticeBolztmann::Muestre(char *nombreArchivo, int t){
  ofstream MiArchivo(nombreArchivo);
  int ix,iy;
  for(ix = 0; ix < Lx ; ix+= 1  ) {
    for( iy = 0 ; iy < Ly ; iy+= 1)
      //cout << ix << " " << iy << endl;
      MiArchivo << ix << " "<< iy << " " << rho(ix,iy,Lz/2, t) <<endl;

  MiArchivo << endl;
  }
  
}
void LatticeBolztmann::Muestre2(int t){
  /* obtener los 
  int ix = 50,iy  ;
  for( iy = 0 ; iy < Ly ; iy++){
    cout << iy << " " << pow(rho(ix , iy , t),2 ) << endl;
  }
  */
  cout << t << " " <<rho(0,0,0,t)<< endl;
}
const int tmax=100;

int main(){
  LatticeBolztmann Ondas;
  
  Ondas.Inicie();
  int t;
  for ( t = 0 ; t < tmax ; t++){
    Ondas.Colisione(t);
    Ondas.Adveccione();
    
        
  }
  
  Ondas.Muestre("D3Q15.dat" , t);
  
  return 0;
}
