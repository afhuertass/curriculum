// Suma de dos vectores c=a+b en CUDA
#include<iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

#define Lx 128
#define Nx 64
const int Mx=(Lx+Nx-1)/Nx;


//--------------- KERNELS ----------------
__constant__ float d_w[5];
__constant__ int d_Vx[5];
__constant__ int d_Vy[5];

__global__ void SumarConstante(float * d_a){
  int ix;
  ix=blockIdx.x*blockDim.x+threadIdx.x;
  d_a[ix]+=d_Vx[1];
}

int main(){
  float h_w[5]; // w[i]
  int h_Vx[5],h_Vy[5]; // Vx[i],Vy[i]
  int ix, iy ;

  //DECLARAR LAS MATRICES
  //CONSTANTES
  //---Cargar las constantes en el Host-----------------
  //Cargar los pesos
  h_w[0]=1.0/3;    h_w[1]=h_w[2]=h_w[3]=h_w[4]=1.0/6;
  //Cargar los vectores
  h_Vx[0]=0;  h_Vx[1]=1;  h_Vx[2]=0;  h_Vx[3]=-1; h_Vx[4]=0;  
  h_Vy[0]=0;  h_Vy[1]=0;  h_Vy[2]=1;  h_Vy[3]=0;  h_Vy[4]=-1;
  //------Enviarlas al Device-----------------
  cudaMemcpyToSymbol(d_w,h_w,5*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vx,h_Vx,5*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vy,h_Vy,5*sizeof(int),0,cudaMemcpyHostToDevice);

  //Datos
  float h_a[Lx];
  float h_b[Lx][Ly];

  //En el Device
  float *d_a;  cudaMalloc((void**) &d_a,Lx*sizeof(float));
  float *d_b; size_t pitcha ; 
  cudaMallocPitch( (void**) &d_a , &pitcha ,Lx*sizeof(float), Lx);

  //INICIALIZAR LOS DATOS
  //Cargar los datos en el Host
  for(ix=0;ix<Lx;ix++)
    h_a[ix]=ix;
  for( ix = 0; ix < Lx; ix++)
    for( iy = 0; iy < Lx ; iy++)
      h_b[ix][iy] = Lx*ix+iy;

  //Enviarlos al Device
  cudaMemcpy(d_b, pitcha , h_b , Lx*sizeof(float) , Lx*sizeof(float) , Lx , cudaMemcpyHostToDevice );
  cudaMemcpy2D(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);
  
  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,1,1);
  dim3 BlocksPerGrid(Mx,1,1);
  SumarConstante<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a);

  //IMPRIMIR LOS DATOS
  //Devolverlos al Host
  cudaMemcpy(h_a,d_a,Lx*sizeof(float),cudaMemcpyDeviceToHost);
  //Imprimirlos
  for(ix=0;ix<Lx;ix++)
    cout<<ix<<" "<<h_a[ix]<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_a);

  return 0;
}
