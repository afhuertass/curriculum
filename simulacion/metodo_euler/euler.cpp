#include <iostream>
#include <cmath>
using namespace std;
double f(double x, double y){
  return y;
}
void UnPasoDeEuler(double &x,double &y, double h){
  //calcular k;
  double k = f(x,y)*h;
  x+=h;
  y+= k;
 
  
}
int main(){
  double x = 0, y = 1;
  double h = 0.1;
  while( x< 10){
    cout << x << " " << y << endl;
    UnPasoDeEuler(x,y,h);
      
  }
}
