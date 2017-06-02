// runge kutta
#include <iostream>
#include <cmath>
using namespace std;
const double omega = 1.1;
double f2(double x, double y1 , double y2){
  return -omega*omega*y1; 
}

double f1(double x, double y1, double y2){
  // funcion a resolver
  return y2;
  
}
double PasoRungeKutta(double &x  , double &y1, double &y2,double h ){
  // definimos cuatro variables k
  double k11,k12, k13,k14;
  double k21,k22, k23,k24;
  double inc = 0, inc2 = 0;
  k11 = h*f1(x,y1,y2);                             k21 = h*f2(x,y1,y2);
  k12 = h*f1(x + h/2 , y1 + k11/2, y2 + k21/2);    k22 = h*f2(x + h/2 , y1 + k11/2 , y2 + k21/2);  
  k13 = h*f1(x + h/2 , y1 + k12/2 , y2 + k22/2);   k23 = h*f2(x + h/2 , y1 + k12/2 , y2 + k22/2 );
  k14 = h*f1(x +h , y1 + k13 ,y2+ k23);            k24 = h*f2(x +h , y1 + k13 ,y2 +k23 );
  
  x+=h;
  //y1+= inc;
  //y2+=inc2;
  y1 += (k11 + 2*k12 + 2*k13 + k14)/6;
  y2 += (k21 + 2*k22 + 2*k23 + k24)/6;  
}
int main(){
  double x = 0, y1 = 0 , y2 = 1;
  double h = 0.1;
  while( x < 10){
    cout <<  x << " " <<  y1 << endl; //PasoRungeKutta(x,y,h);
    PasoRungeKutta(x,y1, y2,h);
    
  }
}
