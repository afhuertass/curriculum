#include<iostream>
#include<cmath>
#include "Vector.h"
// este programa resuelve los puntos b,c,d de la tarea
/* agregando una dos funciones a la calse Cuerpo. 
la primera llamada DibujeseRelativo, recibe el sol y jupiter.
con el sol en el origen, calcula el angulo entre jupiter y el sol en el sistema sin rotar, y con dicho angulo calcula las componentes del nuevo vector de posicion (Del Cuerpo), como aplicando una matriz de rotacion y procede a imprimir dichas coordenadas de forma de gnuplot pueda dibujarlas 

La segunda funcion, ImprimirRotado realiza el mismo procedimiento, pero imprime las coordenadas x o y, esto con el fin de obtener la grafica X Y de las coordenadas del troyano. Ligeras modificaciones de esta funcion permiten obtener o bien el par x y, o bien la coordenada x para obtener la grafica de t vs X 
*/
using namespace std;

const double Deltat=10;
const double G=1;
const double chi=0.193183325037836;
const double Um2chi=1-2*chi;

const int N=3;
//-------------------------------------------------------------
class Cuerpo;
class Colisionador;
//-------------------------------------------------------------

class Cuerpo{
private:
  vector3D r,V,F;  double m,R; //Posicion, Vel, Frueza, masa y Radio

public:
  void Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0);
  void BorreFuerza(void);
  void IncrementeFuerza(vector3D F0);
  void Mueva_r1(double dt);
  void Mueva_V(double dt);
  void Mueva_r2(double dt);
  void Dibujese(void);
  void DibujeseRelativo( Cuerpo sol, Cuerpo Jupiter );
  double ImprimirRotado(Cuerpo sol, Cuerpo Jupiter);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  vector3D GetR(void){ return r; };
  friend class Colisionador;
};   // <-------- TIENE QUE COLOCARSE ESE ;!!!!!!!!!!

void Cuerpo::Inicie(double x0, double y0, double Vx0, double Vy0, double m0, double R0){ //No es una funcion suelta sino que pertence a Cuerpo
  r.cargue(x0,y0,0);  V.cargue(Vx0,Vy0,0);
  m=m0, R=R0;
}
void Cuerpo::BorreFuerza(void){
  F.cargue(0,0,0);
}
void Cuerpo::IncrementeFuerza(vector3D F0){
  F+=F0;
}
void Cuerpo::Mueva_r1(double dt){
  r+=V*(chi*dt); 
}
void Cuerpo::Mueva_V(double dt){
  V+=F*(dt/(2*m));
}
void Cuerpo::Mueva_r2(double dt){
  r+=V*(Um2chi*dt);
}
void Cuerpo::Dibujese(void){
  cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
void Cuerpo::DibujeseRelativo(Cuerpo sol, Cuerpo jupiter){
  // esta funcion recibe el sol y a jupiter, calcula el angulo entre los dos
  // y rota el vector posicion del presente cuerpo en dicho angulo
  // de esta forma el problema es mas general y podrian agregarse mas planetas 
  // y observar su movimiento relativo al sistema sol-jupiter
  /* jupiter, Sol, troyano.
     
     Jupiter.DibujeseRelativo ( sol  , jupiter);
     Sol.DibujeseRelativo ( sol, jupiter);
     
     troyano.DibujeseRelativo (sol , jupiter); 
   */ 
  vector3D dif = jupiter.GetR() - sol.GetR();
  double dx = jupiter.Getx() - sol.Getx();
  double dy = jupiter.Gety() - sol.Gety();
  double tan = dy/dx;
  double angle;
  if( dx == 0 && dy == 0){
    cout << ", " <<0 << "+" <<R<<"*cos(t) ,"<<0<<"+"<<R<<"*sin(t)";  
    return ;
  }
  if (dx < 0 ){ // una correccion 
    tan = tan*-1;
    angle = atan(tan);
    angle = M_PI - angle;
      
  }else {
    angle = atan(tan);
  }
  vector3D dif2 = r - sol.GetR();
  /* el planeta1 es el sol y el vector dif2 es un vector que une al sol con el
     planeta en cuestion, lo siguiente es obtener las coordenadas en el sistema rotante de dicho vector. 
   */
  double x = dif2.x()*cos(angle) + dif2.y()*sin(angle);
  double y = -dif2.x()*sin(angle) + dif2.y()*cos(angle);
  //
  cout << ", " <<x << "+" <<R<<"*cos(t) ,"<<y<<"+"<<R<<"*sin(t)";        
}
double Cuerpo::ImprimirRotado(Cuerpo planeta1 , Cuerpo planeta2){
  // Funcion identica a la anterior pero con el objetivo de retornar la coordenada x en el sistema rotante, o imprimir las coordenas x y en el sistema rotante.
  vector3D dif = planeta2.GetR() - planeta1.GetR();
  double dx = planeta2.Getx() - planeta1.Getx();
  double dy = planeta2.Gety() - planeta1.Gety();
  double tan = dy/dx;
  double angle;
  if( dx == 0 && dy == 0){
    //cout << ", " <<0 << "+" <<R<<"*cos(t) ,"<<0<<"+"<<R<<"*sin(t)";  
    cout << "0 0" << endl;
    return 0.0;
  }
  if (dx < 0 ){
    tan = tan*-1;
    angle = atan(tan);
    angle = M_PI - angle;
      
  }else {
    angle = atan(tan);
  }
  vector3D dif2 = r - planeta1.GetR();
  double x = dif2.x()*cos(angle) + dif2.y()*sin(angle);
  double y = -dif2.x()*sin(angle) + dif2.y()*cos(angle);
  //
  //cout << ", " <<x << "+" <<R<<"*cos(t) ,"<<y<<"+"<<R<<"*sin(t)";
  // cout << x << " " << y << endl; 
  return x;
}

//------------------------------------------------------------------
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Planeta);
  void AgregueFuerza(Cuerpo & Planeta1, Cuerpo & Planeta2);
};

void Colisionador::CalculeFuerzas(Cuerpo * Planeta){
  int i,j;
  for(i=0;i<N;i++)
    Planeta[i].BorreFuerza();
  for(i=0;i<N;i++)
    for(j=0;j<i;j++)
      AgregueFuerza(Planeta[i],Planeta[j]);
}
void Colisionador::AgregueFuerza(Cuerpo & Planeta1, Cuerpo & Planeta2){
  vector3D dr=Planeta1.r-Planeta2.r;
  double aux=G*Planeta1.m*Planeta2.m*pow(norma2(dr),-1.5);
  vector3D F2=dr*aux;
  Planeta2.IncrementeFuerza(F2);  Planeta1.IncrementeFuerza(F2*(-1));
}

//--------------------FUNCIONES GLOBALES-----------------------
void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl;
  //cout<<"set output 'animacion_troyano.gif'"<<endl;
  cout<<""<<endl;
  cout<<"unset key"<<endl;              // QUITAR LABEL
  cout<<"set xrange [-120:120]"<<endl;     // RANGO EN X
  cout<<"set yrange [-120:120]"<<endl;  
  cout<<"set size ratio -1"<<endl;      // VENTANA CUADRADA
  cout<<"set parametric"<<endl;         // GRAFICA PARAMETRICA
  cout<<"set trange [0:7]"<<endl;       // RANGO 
  cout<<"set isosamples 12"<<endl;      // EL BALON ESTA HECHO DE 12 LINEAS

}
void InicieCuadro(void){
  cout<<"plot 0,0 ";
}
void TermineCuadro(void){
  cout<<endl;
}



int main(){
  Cuerpo Planeta[N];
  Colisionador Newton;
  double t,tdibujo;
  double r=100, m0=1047, m1=1;  // datos de entrada para calcular la condcicion inicial
  double r0, r1, V0,V1;
  double M,omega,T;
  int i;
  double mt = 0.005;
  double ang = 60*M_PI/180.0;
  double ang2 = 30*M_PI/180.0;
  M=m0+m1;   omega=sqrt(G*M*pow(r,-3));  T=2*M_PI/omega;  
  r1=m0*r/M;  r0=r1-r;   V1=omega*r1;   V0=omega*r0;

  InicieAnimacion(); // comentamos el codigo referente a la animacion

  Planeta[0].Inicie(r0, 0 ,0 ,V0 ,m0 , r/10); // sol
  Planeta[1].Inicie(r1, 0 ,0 ,V1 ,m1 , r/20); // jupiter
  Planeta[2].Inicie(r1*cos(ang)+0.01, r1*sin(ang)+0.01 ,-V1*cos(ang2),V1*sin(ang2),mt , r/30 ); 

// el troyano, con la misma velocidad de jupiter y corrido 60 grados respecto al sol
  for(t=tdibujo=0;t<40*T;t+=Deltat,tdibujo+=Deltat){
    
    /* if(tdibujo>T/500){
      
    
     Planeta[1].DibujeseRelativo(Planeta[0] , Planeta[1]); // dibuje jupiter relativo al sol 
    Planeta[0].DibujeseRelativo(Planeta[0] , Planeta[1]);
    Planeta[2].DibujeseRelativo(Planeta[0], Planeta[1]); // troyano. 
    
        tdibujo=0;
    }*/
  // Calculamos con el metodo de verlet optimizado
     InicieCuadro();
     for ( i = 0 ; i< N ; i++) {
       Planeta[i].DibujeseRelativo( Planeta[0], Planeta[1]);
       Planeta[i].Dibujese();
     }

    TermineCuadro();


    for(i=0;i<N;i++) Planeta[i].Mueva_r1(Deltat);
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<N;i++){Planeta[i].Mueva_V(Deltat);  Planeta[i].Mueva_r2(Deltat);}
    Newton.CalculeFuerzas(Planeta);
    for(i=0;i<N;i++){Planeta[i].Mueva_V(Deltat);  Planeta[i].Mueva_r1(Deltat);}

    /* Descomentando las lineas de la animacion es posible obtener un gif como el que se adjunta, con los dos sistema mostrados, el rotado y sin rotar  
    
       

     */ // para obtener la posicion X del troyano en funcion del tiempo se descomentan las siguientes dos lineas 
	 
    //double xr = Planeta[2].ImprimirRotado( Planeta[0] , Planeta[1]);
    // cout << t << " " << xr << endl;
  }

  return 0;
}


