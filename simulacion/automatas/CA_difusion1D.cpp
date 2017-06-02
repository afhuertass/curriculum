#include <iostream>
#include <cmath>

#include "Random.h"

using namespace std;

const int Lx = 100;
const double p = 0.5;


class Automata{
private:
  int **n; 
  int **nnew;
  int Vx[2];
public: 
  Automata(void);
  ~Automata(void);
  void Inicie(double mu, double sigma, int N, Crandom & ran2);
  void Evolucione(Crandom & ran2);
  int Getn(int ix){ return n[ix][0] + n[ix][1]; };
  void Muestre(void);
};
Automata::Automata(void){
  n = new int*[Lx]; nnew = new int*[Lx]; // apuntadores
  for (int i = 0 ; i< Lx; i++){
    n[i] = new int[2]; 
    nnew[i] = new int[2]; 
    
  } // con
  Vx[0] = 1; Vx[1] = -1;
  
}
void Automata::Inicie(double mu, double sigma, int N, Crandom & ran2){
  // borrar todos los contenidos
  int ix, j , c;
  for ( ix = 0; ix < Lx ;ix++){
    for (j =0 ;j<2;j++){
      n[ix][j] = 0;
    }
  }
  // 
  c = 0;
  while(c < N/2){
    ix = (int) ran2.gauss(mu,sigma);
    if( ix >=0 && ix < Lx){
      if(n[ix][0] == 0) n[ix][0] = 1 ; c++; 
      
    }
  }

  // a la izquierda
  c = 0;
  while(c < N/2){
    ix = (int) ran2.gauss(mu,sigma);
    if( ix >=0 && ix < Lx){
      if(n[ix][0] == 0) n[ix][1] = 1 ; c++; 
      
    }
  }
}
void Automata::Evolucione(Crandom & ran2){
  int ix,j;
  int **aux;
  double r;
  for(ix =0, ix < Lx; ix++){ // colision
    r = ran2.r();
    if(r < p){
      // quede igual
      for(j = 0; j< 2 ; j++){
	nnew[ (ix + Vx[j]+Lx)%Lx ][j] = n[ix][j];
      }
    }else { // interecambiar
       for(j = 0; j< 2 ; j++){
	 nnew[ (ix + Vx[j]+Lx)%Lx ][j] = n[ix][j+1]%2;
      }
      
    }
  }
  
  aux = n;
  n = nnew;
  nnew = aux;

  
}
void Automata::Muestre(void){
  int j,ix;
  for (j = 0 ; i<2 ; j++){
    for (ix = 0 ; ix < N ; ix++){
      cout << n[ix][j];
      cout << endl;
    }
  } 
}
Automata::~Automata(void){
  for (int i = 0 ; i< Lx; i++){
    delete[] n[i];
    delete[] nnew[i];
    
  } // con
  delete[] n; delete[] nnew[]
}


int mai(void){
  Automata difusion;
  Crandom ran2(10);
  
  // difusion.Inicie( 0.5*Lx, 0.1*Lx, 2 ,ran2);
  

  return 0;
}
