#include <iostream>
#include <cmath>
// implentacion modelo sir. ecuaciones acopladas. 
using namespace std;
// parametros del modelo 



double f1(double r,double r1, double r2){ // funcion ds/dt = f1 = -betasi
  return r2;
}
double f2(double r, double r1, double r2, double lambda){
  double lam2 = lambda*lambda;
  double r_2 = r*r;
  return (r1*lam2/r_2 -r1 -r2/r);

}

double PasoRungeKutta(double &r , double &r1 , double &r2 , double h, double lam){
double k11,k12, k13,k14;
double k21,k22, k23,k24;

 k11 = h*f1(r,r1,r2);                   k21 = h*f2(r,r1,r2,lam);
 k12 = h*f1(r+h/2 , r1+k11/2, r2+k21/2); k22 = h*f2(r+h/2 , r1+k11/2,r2+k21/2,lam);
 k13 = h*f1(r+h/2 , r1+k12/2 ,r2+k22/2); k23 = h*f2(r+h/2 , r1+k12/2,r2+k22/2,lam);
 k14 = h*f1(r+h, r1+k13, r2+k23);         k24 = h*f2(r+h, r1+k13, r2+k23,lam);

r+=h;
r1+= (k11 + 2*k12 + 2*k13 + k14)/6;
r2+= (k21 + 2*k22 + 2*k23 + k24)/6;

}
int main(void){
  double r = 0.01,  r1 = 1 , r2 = 0;
  double h = 0.1;
  while ( r < 10){
   
// tiempo suceptibles infectados retirados
    cout << r << " " << r1 << endl;
    PasoRungeKutta(r,r1,r2,h,0);
  }
}
