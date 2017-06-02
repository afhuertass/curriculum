#include <iostream>
const int Lx = 100, Ly=100 , Lz = 100;
const int q15 = 15;
using namespace std;
int main()
{
  int i = 0,j = 0 , k = 0, q = 0;
  double**** f;
  double**** fnew;
  
  f = new double***[Lx];
  fnew = new double***[Lx];
  for(i = 0; i < Lx ; i++){
    f[i] = new double**[Ly];
    fnew[i] = new double**[Ly];
    for (j = 0 ; j < Ly ; j++){
      f[i][j] = new double* [Lz];
      fnew[i][j] = new double*[Lz];
      for ( k = 0 ; k < Lz ; k++){
	f[i][j][k] = new double[q15] ;
	fnew[i][j][k] = new double[q15];
      }
	
    }
  }
  f[0][0][0][0] = 10 ;
  cout << f[0][0][0][0] << endl;
}
