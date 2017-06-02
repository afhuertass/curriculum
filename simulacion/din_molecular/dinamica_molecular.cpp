
#include <iostream>
#include <cmath>

using namespace std;
const double g=9.8 //  gravedad m/ss
  const double Deltat = 0.1;
  class Cuerpo{
  private:
    double x,y,Vx,Vy, Fx,Fy,m,R;
  public:
   void Inicie(double x0,double y0, double Vx0, double Vy0 , double m0, double R0);
    void CalcularFuerza(void);
    void Muevase(double dt);
    
  };
void Cuerpo::Inicie(double x0,double y0, double Vx0, double Vy0 , double m0, double R0){
  x = x0; y = y0; 
  Vx = Vx0; Vy0 = Vy0;
  m = m0;
  R = R0;
}
void Cuerpo::CalcularFuerza(void){
  Fx=0; Fy=-m*g;

}
void Cuerpo::Muevase(double dt){
  x += Vx*dt;y += Vy*dt;
  Vx += Fx*dt/m ; Vy += Fy*dt/m;
}
int main(){
  Cuerpo balon;
  balon.Inicie(0,0,3,4,0.453,0.15);
  for (double t = 0; t < 10 ; t += Deltat){
    cout << balon.x << " " balon.y << endl;
    balon.CalculeFuerza();
    balon.Muevase(Deltat);
    
  }
  return 0;
  
}
