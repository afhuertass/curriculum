#include <iostream>

int  main(){

  double under = 1;
  double over = 1;
  int n = 100;
    std:: cin >> n;
  for(int i=1 ; i <n; i++){
    under = under/2.0;
    over = over*2.0;
    std::cout << "Under " << under << "\n" ;
    std::cout << "Over " << over << "\n";

  }
  return 0;

}
