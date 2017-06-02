//Mi Primer Programa
#include <iostream>
#include <cmath>
using namespace std;
const double ERR = 1e-7;
double f(double x){
  return sin(x)/x;

}
int main(){
  //Defino a y b
  double a = 2 , b = 4;
  double fa = f(a), fb = f(b) , fm; // calcular fa y fb
  double m;
  while( (b-a)> ERR){
    m = (a+b)/2;
    fm = f(m);
    if( fa*fm < 0 ){ // el corte esta en la primera parte
      b = m ;
      fb = fm;
    }else {
      a = m;
      fa = fm;
    }
  }
  cout << "el cero esta en: " << (a+b)/2 << endl;
  cout << "El error fue:" << (a+b)/2-M_PI << endl;
}
