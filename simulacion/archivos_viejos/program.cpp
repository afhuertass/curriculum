#include <iostream>

int main(int arg, char* arg2[]){
  
  for(int i = 0; i<10; i++){
    
    for(int j = 0; j < i; j++){
      std::cout << "*";

    }
    std::cout << "\n";
  }
  return 0;
}
