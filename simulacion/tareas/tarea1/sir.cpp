#include <iostream>
#include <cmath>
// implentacion modelo sir. ecuaciones acopladas. 
using namespace std;
// parametros del modelo 
const double Pbeta = 0.5;
const double Pgamma = 0.08;

double f1(double t,double s, double i){ // funcion ds/dt = f1 = -betasi
return -1*Pbeta*s*i;
}
double f2(double t, double s, double i){ // funcion di/dt = b si-gamma i

return Pbeta*s*i-Pgamma*i;

}
double Retirados(double s, double i){
// s + i + r = 1 
return (1-s-i);

}
double PasoRungeKutta(double &t , double &s , double &i , double h){
double k11,k12, k13,k14;
double k21,k22, k23,k24;

k11 = h*f1(t,s,i);                     k21 = h*f2(t,s,i);
k12 = h*f1(t+h/2 , s+k11/2, i+k21/2);  k22 = h*f2(t+h/2 , s+k11/2, i+k21/2);
k13 = h*f1(t+h/2 , s+k12/2 , i+k22/2); k23 = h*f2(t+h/2 , s+k12/2 , i+k22/2);
k14 = h*f1(t+h, s+k13, i+k23);         k24 = h*f2(t+h, s+k13, i+k23);

t+=h;
s+= (k11 + 2*k12 + 2*k13 + k14)/6;
i+= (k21 + 2*k22 + 2*k23 + k24)/6;

}
int main(void){
double t = 0,  s = 0.999 , i = 0.001 , r = 0;
double h = 1;
while ( t < 100){
  r= Retirados(s,i); 
  // tiempo suceptibles infectados retirados
  cout << t <<  " " << s << " " << i << " " << r << endl;
  PasoRungeKutta(t,s,i,h);
}
}
