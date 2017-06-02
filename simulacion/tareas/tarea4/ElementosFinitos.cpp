
#include<iostream>
#include<cmath>
#include"LU.h"

using namespace std;

const int A=1.0; //[0:A]
const int E=4.0; //Valor minimo 3.

const int M=E-1.0; //# Elementos
const int S=E-2.0;
const double dX=(double) A/M;

//Condiciones de Frontera
const double dUdx0=1.0;
const double dUdxA=1.0;

class Elemento{
private:
  double  Epsilon[E][E][E];//Matriz[Elemento][Fila][Columna]
  double  RHS[E],U[E],UF[S],RHSF[S];
  double  **EpsilonTotal;
  double  **EpsilonFinal;
  
public:
  Elemento(void);
  ~Elemento(void);
  void InicieElemento(void);
  void CalculeEpsilon(void);
  void ResuelvaSistema(void); 
};

//------------------------------------------

Elemento::Elemento(void){
  EpsilonTotal=new double* [E];
  for(int ix=0;ix<E;ix++)
    {EpsilonTotal[ix]=new double [E];}

  EpsilonFinal=new double* [S];
  for(int ix=0;ix<S;ix++)
    {EpsilonFinal[ix]=new double [S];}
}


Elemento::~Elemento(void){
  for(int ix=0;ix<E;ix++)
    {delete[] EpsilonTotal[ix];}
  delete[] EpsilonTotal;

 for(int ix=0;ix<S;ix++)
    {delete[] EpsilonFinal[ix];}
  delete[] EpsilonFinal;
}

void Elemento::InicieElemento(void){
  int i,j,k;

 for(k=0;k<E;k++)
  for(i=0;i<E;i++)
    for(j=0;j<E;j++)
      Epsilon[k][i][j]=0;

 for(i=0;i<E;i++)
    for(j=0;j<E;j++)
      EpsilonTotal[i][j]=0;

 for(i=0;i<E;i++){
   RHS[i]=0;
  if(i==0)  RHS[i]=-dUdx0;
  if(i==E-1)RHS[i]=dUdxA; 
  }
 //Condiciones de Frontera:
 U[0]=0; U[E-1]=A;
}

void Elemento::CalculeEpsilon(void){
  int i,j,k;

  for(k=0;k<E;k++)
  for(i=0;i<E;i++)
    for(j=0;j<E;j++)
      if(i==k && j==k)     Epsilon[k][i][j]=(1.0/dX)+(dX/3.0);
      else if(i==k && j==k+1)   Epsilon[k][i][j]=(-1.0/dX+dX/6.0);
      else if(j==k && i==k+1)   Epsilon[k][i][j]=(-1.0/dX+dX/6.0);
      else if(i==k+1 && j==k+1) Epsilon[k][i][j]=(1.0/dX)+(dX/3.0);
   

 for(k=0;k<E-1;k++) 
  for(i=0;i<E;i++)
    for(j=0;j<E;j++)
      EpsilonTotal[i][j]+=Epsilon[k][i][j];  

 for(i=0;i<S;i++) 
   for(j=0;j<S;j++)
      EpsilonFinal[i][j]=EpsilonTotal[i+1][j+1];
      
 for(i=0;i<S;i++) 
   {   RHSF[i]=RHS[i+1]; if(i==S-1){RHSF[i]=-EpsilonTotal[i][i+1]; } }
 // RHSF[i]=-EpsilonTotal[1][i+1];
 

 /*
 //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
   cout<<"ET"<<endl;      

  for(i=0;i<E;i++){
	cout<<"["; 
      for(j=0;j<E;j++)
	{cout<<EpsilonTotal[i][j]<<" ";}
  	cout<<"]"; 
       cout<<endl; }
 cout<<endl;      
 cout<<endl;      
   //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 cout<<"EFinal "<<endl;      
   for(i=0;i<S;i++){
     cout<<"["; 
      for(j=0;j<S;j++)
	{cout<<EpsilonFinal[i][j]<<" ";}
  	cout<<"]"; 
       cout<<endl; }
   cout<<endl;  
   cout<<endl;          
   //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  cout<<"RHS"<<endl;
  for(i=0;i<S;i++) 
     cout<<RHSF[i]<<" "<<EpsilonTotal[1][i+1]<<endl;

 */
  
}

void Elemento::ResuelvaSistema(void){
  //LU: Ax=b: Epsilon U = RHS.

  ResolverLU(EpsilonFinal,RHSF,UF,S);
 int i,j;
 for(i=0;i<S;i++)
   if(i<E-1)
     U[i+1]=UF[i]; 

 //Imprimir 
 // cout<<"[Epsilon]"<<" "<<"[U]"<<"="<<"[RHS]"<<endl;
  cout<<endl;
  for(i=0;i<S;i++){
	cout<<"["; 
      for(j=0;j<S;j++)
	{cout<<EpsilonFinal[i][j]<<" ";}
	cout<<"] "<<UF[i]<<"="<<RHSF[i];
       cout<<endl; }
  cout<<endl;
 
 /*
  cout<<endl;  
   cout<<endl; 
  for(i=0;i<E;i++){
	cout<<"["; 
      for(j=0;j<E;j++)
	{cout<<EpsilonTotal[i][j]<<" ";}
	cout<<"] "<<U[i]<<"="<<RHS[i];
       cout<<endl; }
   cout<<endl;
   cout<<endl;  
 */
  double R;
 
 cout<<"Posicion"<<" "<<"Aprox"<<" "<<"Exacta"<<" "<<"Error"<<endl;
  for(i=0;i<S;i++){  
    R=sinh((i+1)*dX)/sinh(1.0);
    cout<<(i+1)*dX<<" "<<UF[i]<<" "<<R<<" "<<(R-UF[i])/R*100<<endl;
  
    }
 
}



 int main(void){

  Elemento Finito;

  Finito.InicieElemento();
  Finito.CalculeEpsilon();
  Finito.ResuelvaSistema();

  return 0;
}
