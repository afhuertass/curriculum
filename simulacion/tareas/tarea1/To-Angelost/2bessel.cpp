#include <iostream>
#include <cmath>
// Aqui se soluciona la ecuacion de bessel utilizando un runge-kutta de 4 orden
// donde se puede obtener la solucion para cualquier lambda, con una ligera modificacion se puede conseguir la solucion para cualquier lambda, en particular modificando el main para que no itere sobre muchos valores de lambda, se obtuvo la solucion para la parte a.

// Este programa tambien calcula las raices de f(lambda)    
using namespace std;

double f1(double r,double r1, double r2){ // funcion ds/dt = f1 = -betasi
  return r2;
}
double f2(double r, double r1, double r2, double lambda){
  // la solucion depende el parametro lambda 
  double lam2 = lambda*lambda;
  return (-lam2*r1-r2/r );

}

double PasoRungeKutta(double  &R , double &R1 , double &R2 , double H, double LAM){
double K11,K12, K13,K14;
double K21,K22, K23,K24;

 K11 = H*f1(R,R1,R2);                   K21 = H*f2(R,R1,R2,LAM);
 K12 = H*f1(R+H/2 , R1+K11/2, R2+K21/2); K22 = H*f2(R+H/2 , R1+K11/2,R2+K21/2,LAM);
 K13 = H*f1(R+H/2 , R1+K12/2 ,R2+K22/2); K23 = H*f2(R+H/2 , R1+K12/2,R2+K22/2,LAM);
 K14 = H*f1(R+H, R1+K13, R2+K23);         K24 = H*f2(R+H, R1+K13, R2+K23,LAM);

R+=H;
R1+= (K11 + 2*K12 + 2*K13 + K14)/6;
R2+= (K21 + 2*K22 + 2*K23 + K24)/6;
 return 0.0;
}
int main(void){
  double r = 0.01,  r1 = 1 , r2 = 0 ;
  double h = 0.1, lambda = 0.01; // empezamos con lambda en 0.01 
  int size = 15/0.001, i = 0 ; //  aumentando en pasos de 0.001 hasta 15, para
  //tener muchos puntos
  double lambdas[size], flam[size];
  // en estos arreglos guardamos los lambdas y el valor de la funcion en r = 1
  
  while ( lambda <= 15){ // iteramos sobre lambda
    r = 0.01 ; r1 = 1 ; r2 = 0; // definimos las condiciones iniciales sobre lambda
    while ( r < 10){
      // este bucle aplica el metodo runge-kutta 
      if( r == 1.01){  // en r = 1 
	//cout << lambda << " " << r1 << endl; 
	lambdas[i] = lambda; // guardamos las parejas lambda y f(lambda)
	flam[i] = r1;
	i++;
      }
      cout << r << " " << r1 << endl;
      PasoRungeKutta(r,r1,r2,h,lambda);
    }
    lambda +=100;
  }
  // Calculamos las raices 
  int roots = 0; // esta variable es para enumerar las raices encontradas. 
  for(int j = 0; j < (size-1) ; j++ ){
    double lam1, lam2;
    double f1, f2;
    lam1 = lambdas[j]; //obtenemos f en la entrada j  
    lam2 = lambdas[j+1]; // y f en la entrada siguiente j+1
    // <---
    f1 = flam[j];
    f2 = flam[j+1];
    if( f1*f2 < 0 ){ // si tienen signo opuesto la raiz esta entre los dos
      
      // finalmente estimamos que la raiz esta en el promedio de los dos lambdas
      cout <<"#"<<  roots << " " << (lam1+lam2)*0.5 << endl;
      
      roots++;
    }
    
  }
}
