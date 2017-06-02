// runge kutta
#include <iostream>
#include <cmath>
using namespace std;
double f(double x, double y){
  // funcion a resolver
  return y;
  
}
double PasoRungeKutta(double &x  , double &y, double h){
  // definimos cuatro variables k
  double k1,k2, k3,k4;
  double inc = 0;
  k1 = f(x,y);
  k2 = f(x + h/2 , y + k1*h/2);
  k3 = f(x + h/2 , y + h*k2/2);
  k4 = f(x +h , y + h*k3);
  //cout << k1 << " " << k2 << " " << k3 << endl;
  inc = k1+2*k2+2*k3 + k4;
  inc = h/6*inc;
  //cout << inc << endl;
  x+=h;
  y+= inc;
}
int main(){
  double x = 0, y = 1;
  double h = 0.1;
  while( x < 10){
    cout <<  x << " " <<  y << endl; //PasoRungeKutta(x,y,h);
    PasoRungeKutta(x,y,h);
    
  }
}
