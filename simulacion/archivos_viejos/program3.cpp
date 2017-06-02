#include <iostream>

int main(int arg, char* args[]){
  int col = 5;
  int interval = 0;
  for(int i = 0; i < 6;i++){
    
    for(int j = 0; j< 11 ; j++){
      if(j >= col-interval && j <= col+interval  ){
	if(i%2 == 0){
	  std::cout << "*";
	}else {
	std::cout << "a";
	}
      }else {
	std::cout << " ";
      }
    }
    interval++;
    std::cout << "\n";

  }
  return 0;

}
