//Mi Primer Programa
#include <iostream>
#include <cmath>
using namespace std;

double f(double alpha, double x,double t){
  return cos(alpha*t-x*sin(t));

}

double IntegralPorSimpson(double alpha,double x , double a,double b, int N){
  double h,suma; int i;
  double t;
  h=(b-a)/(2*N);
  suma=0;
  for(i=0;i<=2*N;i++){
    t=a+i*h;
    if(i==0 || i==2*N)
      suma+=f(alpha,x,t);
    else if(i%2==0)
      suma+=2*f(alpha,x,t);
    else
      suma+=4*f(alpha,x,t);
  }
  return h/3*suma;
}
double Jbessel(double alpha,double x){
  return 1.0/(2*M_PI)*IntegralPorSimpson(alpha,x,0,2*M_PI , 50);
}
int main(){
  double a=0,b=10; int N=5;
  double  alpha = 0 , x ;
  //for (alpha = 1 ; alpha < 10 ; alpha++) {
    for(x =0 ; x< 10 ; x+=0.1){
      cout << x << " " << Jbessel(alpha,x) << endl;
    }
    //}
  //cout<<"La integral vale ="<<IntegralPorSimpson(a,b,N)<<endl;

  return 0;
}
