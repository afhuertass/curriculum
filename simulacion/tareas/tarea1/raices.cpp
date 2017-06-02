#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
using namespace std;
// este programa cargara 
int main(void){
  int n_lines = 0;
  string line;
  ifstream myfile ("2bessel.dat");
  if( myfile.is_open() ){
    
    while( getline(myfile,line)){
      
      cout << line << endl;
    }
    myfile.close();
  }

  

}
