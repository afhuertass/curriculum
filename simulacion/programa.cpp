// coment!
#include <iostream>
#include <cmath>

using namespace std;
double f(double x){
  return x*x;
}

int main(){
 
  double x,h,suma; 
  double a =0, b = 10;
  int N = 5 ;
  int i = 0;
  h = (b-a)/2*N;
  
  for( i = 0 ;  i<2*N; i++){
    //cout << 2*i+1 << endl;
    //suma += 2*i+1;
    //cout << x << " " << pow(x,1.5) << endl;
    x = a + i*h;
    if ( i == 0 || i ==2*N){ // primero o el ultimo
      suma += f(x);
    else if 
      if(i%2==0 ){ // es par
	suma += 2*f(x)
      }else { // es impar
	suma +=4*f(x)
      }
    
  }
  cout << "la integral es:" << suma << endl;
  return 0;
}



