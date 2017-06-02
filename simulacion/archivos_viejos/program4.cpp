#include <iostream>
float factorial(int);
float potencial(float, int);
int main(int arg, char* args[]){
  
 float x = 0.1;
 int n = 20;
  
  float res = 0;
  for(int j = 0; j < 20; j++){
    res += potencial(x,j)/factorial(j);
  }
  std::cout << "e a la 0.10 = " << res << "\n";
}

float factorial(int n){

  float parc = 1;
  if(n == 0){
    return 1;
  }
  for(int i = 1; i <=n ; i++){
    parc *= i; 
  }
  return parc;
}
float potencial(float x, int n ){
  float resu=x;
  
  for(int i = 1 ; i<=n; i++){
    resu *= x;

  }
  return resu;
}
