// Con el fin de mejorar el programa
// se extendera de memoria stack a memoria dinamica, evitando asi el over flow
// las variables f y fnew deberan ser arreglos dinamicos

#include <iostream>
#include <fstream>
#include <cmath>
// #include "Random64.h"
using namespace std;
// Automata difusion continua
const int Lx=100 ,Ly = 150 ; // tamaño cuadricula 
const int Lz =100;
const double C =  0.5;// C < sqrt( 3/5)
const double C2 =C*C;
const double Um3Umw0  = 1-2*C2;
const double tau = 0.5;

const double A = 2000;
const double lambda = 2;
const double omega =(2*M_PI*C)/lambda; 
const int R = 20;
const int lx =  Lx/2, ly = 4 + 5 + R , lz = Lz/2; 

enum tipoCelda { normal,  abso, fuente, espejo, lente};

class LatticeBolztmann{
private:
  //double f[Lx][Ly][Lz][15], fnew[Lx][Ly][Lz][15]; // 
  double**** f; // en el constructor del lattice 
  double**** fnew;
  
  int V[3][15]; 
  double w[15] ; // pesos de la funcion de equilibrio. 
  tipoCelda CeldaTipo[Lx][Ly][Lz];  // cero = normal , 1 = absorbente,  2 = fuente , 3 = espejo 

public:
  LatticeBolztmann(void);
  ~LatticeBolztmann();
  void Inicie(void);
  void Muestre(char *nombreArchivo,int t);
  void Muestre2(int t);
  void Muestre3(int t);
  double rho(int ix, int iy,int iz, int t);
  double Jx(int ix, int iy, int iz);
  double Jy(int ix, int iy ,int iz);
  double Jz(int ix,int iy , int iz);
  double feq(double rho0, double Jx0 , double Jy0 , double Jz0 , int i , double V ); 
  double SacarVelocidad(int ix, int iy, int iz);
  void Colisione(int t);
  void Adveccione(void);
  void ConstruyaGeometria(void);
  };
double LatticeBolztmann::SacarVelocidad(int ix, int iy, int iz){
  if ( CeldaTipo[ix][iy][iz] == normal) {
    return 0.5;
  }
  if (CeldaTipo[ix][iy][iz] == lente){
    return 0.5/2;
  }
  
}



void LatticeBolztmann::ConstruyaGeometria(void){
  int ix, iy, iz;

  //Lente
  for( ix = 0 ; ix < Lx ; ix++)
    for(iy = 0 ; iy < Ly ; iy++)
      for ( iz = 0 ; iz < Lz ; iz++) {
	CeldaTipo[ix][iy][iz] = normal;
      } 

 // Paredes  
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
  

  //Fuente
  iy=1;
  for(ix=0; ix<Lx; ix++){
    for(iz=0; iz<Lz; iz++){
      CeldaTipo[ix][1][iz] = fuente;
    }
  }
 
  //Rejilla
  iy = 4; 
  double xr = Lx/2, anchox = 5, zr= Lz/2, anchoz=15;
  
  for( ix = 0 ; ix < Lx ; ix++)
    for(iz = 0 ; iz < Lz ; iz++){ 
      if ( (ix > xr - anchox/2.0 && ix < xr + anchox/2.0) && ( iz > zr - anchoz/2 && iz < zr + anchoz/2) ) {
        CeldaTipo[ix][iy][iz] = normal;
      }
     else   {
        CeldaTipo[ix][iy][iz] = abso;
      } 
    
    }
  
int Radio = Ly - 4 - 5;
int ixr = Lx/2 , iyr = Ly, izr = Lz/2;
for ( ix = 0 ; ix < Lx ; ix ++)
  for ( iy = 0;iy <Ly ; iy++ )
    for( iz = 0; iz < Lz ; iz++)
       if (  ( pow( ix - ixr,2) + pow(iy - iyr,2) + pow(iz - izr,2) <= pow(Radio,2) ) && (iy < (Ly-Radio+20) ) ) {
        CeldaTipo[ix][iy][iz] = lente

       }

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
  double rho0, Jx0, Jy0 , Jz0, vel;
  for(ix =0 ; ix < Lx ; ix++ )  //para cada celda
    for(iy =0 ; iy < Ly ; iy++) 
      for( iz = 0 ; iz < Lz ; iz++)
  {
    rho0 = rho(ix,iy,iz,t);
    Jx0 = Jx(ix,iy,iz); 
    Jy0 = Jy(ix,iy,iz);
    Jz0 = Jz(ix,iy,iz);
    vel = SacarVelocidad(ix,iy,iz);
    vel = vel*vel;
    // para cada direccion
    for(i =0 ; i < 15 ; i++)
  // evolucion ecuacion de bolztmann BGK
      fnew[ix][iy][iz][i] = f[ix][iy][iz][i] - 1.0/tau*( f[ix][iy][iz][i] - feq(rho0,Jx0,Jy0,Jz0,i , vel) );  
  }
}

double LatticeBolztmann::feq(double rho0, double Jx0 , double Jy0, double Jz0 , int i, double vel ){
  
  if(i == 0)
    //  return 9.0/2*w[0]*rho0*(1-7.0/3*C2);  //rho0*Um3Umw0;
  return 9.0/2*w[0]*rho0*(1-7.0/3*vel);  //rho0*Um3Umw0;
  else 
    return 3*w[i]*(vel*rho0 + V[0][i]*Jx0 + V[1][i]*Jy0 + V[2][i]*Jz0 );
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
  // 
   w[0] = 2.0/9;
  w[1] = w[2] = w[3] = w[4] = w[5] = w[6] =1.0/9;
  w[7] = w[8] = w[9] = w[10] = w[11] = w[12] = w[13] = w[14] = 1.0/72;
 
  // Vectores del modelo 3d q15 
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

  // alojar espacios de memoria para f y f new
  int i,j,k, q15 = 15;
  f = new double***[Lx];
  fnew = new double***[Lx];
  for(i = 0; i < Lx ; i++){
    f[i] = new double**[Ly];
    fnew[i] = new double**[Ly];
    for (j = 0 ; j < Ly ; j++){
      f[i][j] = new double* [Lz];
      fnew[i][j] = new double*[Lz];
      for ( k = 0 ; k < Lz ; k++){
  f[i][j][k] = new double[q15] ;
  fnew[i][j][k] = new double[q15];
      }
  
    }
  } 
}
LatticeBolztmann::~LatticeBolztmann(){
  // liberar la memoria 
  int i,j,k;
  for( i= 0 ; i <Lx; i++){
    for( j = 0; j < Ly ; j++){
      for( k = 0 ; k < Lz ; k++){
  delete[] f[i][j][k];
  delete[] fnew[i][j][k];
      }
      delete[] f[i][j];
      delete[] fnew[i][j];
    }
    delete[] f[i];
    delete[] fnew[i];
  }
  delete[] f;
  delete[] fnew;
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

void LatticeBolztmann::Muestre(char *nombreArchivo,int t){
  ofstream MiArchivo(nombreArchivo);
  int ix, iy = 12, iz ;
  for(ix = 0; ix < Lx ; ix+= 1  ) {
    for( iz = 0 ; iz < Lz ; iz+= 1)
      //cout << ix << " " << iy << endl;
      MiArchivo << ix << " "<< iz << " " << pow( rho(ix,iy,iz, t) ,2 ) <<endl;
    //cout << ix << ", " << pow( rho(ix,iy,Lz/2, t) ,2 ) <<endl;

  MiArchivo << endl;
  }
  
}
void LatticeBolztmann::Muestre3(int t){
  int ix, iz ; // t no es t sino y
  
  for( ix = 0; ix < Lx ; ix++){
      for(iz = 0 ; iz <Lz ; iz++)
	cout << ix << " " << iz << " " << pow( rho(ix,t,iz, t) ,2 ) << endl;
  
      cout << endl;
  }
}
void LatticeBolztmann::Muestre2(int t){
  ofstream MiArchivo("Lentes.dat");
  int ix = 0,iy= 4,  iz = 0  ;
  for ( ix = 0; ix < Lx ; ix++) {
  for( iz = 0 ; iz < Lz ; iz++){
    int val = 0;
    // cout << ix << " " << pow(rho(ix , iy , iz , t) ,2 ) << endl;
    if ( CeldaTipo[ix][iy][iz] == abso || CeldaTipo[ix][iy][iz] == lente ) val = 1;
    MiArchivo << ix << " " << iz << " " <<  val << endl;
  }

  MiArchivo << endl;
  } 
}

void InicieAnimacion(void){
  cout << "" << endl;
  cout << "set term png" << endl;
  //cout << "set output 'aqui.png'" << endl;
  cout << "set pm3d map" << endl;
  //cout << "set term jpeg" << endl;
  cout << "set xlabel 'X' " << endl;
  cout << "set ylabel 'Z' " << endl;
  
}
void IniciarCuadro(int iy){
  double F = pow(15,2)/(iy*lambda);
  cout << "set title 'F=" << F <<"'"<< endl;
  cout << "set output './images/anima"<< iy << ".png'" << endl;
  cout << "splot '-' " << endl ; 
}
void TerminarCuadro(void){
  cout << "e" <<endl;
  

}
const int tmax=200;
int main(){
  LatticeBolztmann Ondas;
  
  Ondas.Inicie();
  int t;
  for ( t = 0 ; t < tmax ; t++){
    Ondas.Colisione(t);
    Ondas.Adveccione();
    
        
  }
  
  Ondas.Muestre( "Fraun.dat", t);
  //Ondas.Muestre2(t); 
  InicieAnimacion();
  int iy;
  for( iy = 4; iy < Ly ; iy++) {
    IniciarCuadro(iy);
    Ondas.Muestre3(iy);
    TerminarCuadro();
  
  }
  return 0;
}
