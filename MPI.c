#include<stdio.h>
#include<stdlib.h>  
#include <sys/time.h>
#include <mpi.h>
#include <time.h>

void ProcesoTipoA(int, int, int);
void ProcesoTipoB(int, int, int);
void maximoD(double *, int, int, double*);
void minimoA(double *,int, int, double*);
void numeroXMatriz(int, int, int, double *, double *);
void promedio(double *, int, int, double*);
void matmulblks(double *, double *, double *, int, int, int);
void matmulblks1(double *, double *, double *, int, int, int);

double dwalltime(){
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

void ProcesoTipoA(int n, int bs, int cantProcesos){
  double *A, *B, *C, *D, *R, *P, *result1, *result2, *resultABC, *resultDCB;
  int i, j;
  int filas = n/cantProcesos;
  double min, max, suma, tMIN, tMAX, tSUMA;

  // Alocar  
  A = (double *) malloc(n*n*sizeof(double));
  B = (double *) malloc(n*n*sizeof(double));
  C = (double *) malloc(n*n*sizeof(double));
  D = (double *) malloc(n*n*sizeof(double));
  R = (double *) malloc(n*n*sizeof(double));
  P = (double *) malloc(filas*n*sizeof(double));
  result1 = (double *) malloc(filas*n*sizeof(double));
  result2 = (double *) malloc(filas*n*sizeof(double));
  resultABC = (double *) malloc(filas*n*sizeof(double));
  resultDCB = (double *) malloc(filas*n*sizeof(double));  

  // Inicializacion
  
  srand(time(NULL));
  int AN = rand() %100;
  int DN = rand() %100;
  
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
	A[i*n + j] = AN;
	B[j*n + i] = 0;
	C[j*n + i] = 0;
	if(j == i){
	  B[j*n+i] = 1;
	  C[j*n+i] = 1;
	}
	D[i*n + j] = DN;
    }
  }

  for (i = 0; i < filas; i++){
    for (j = 0; j < n; j++){
	P[i*n + j] = 0;
	R[i*n + j] = 0;
	result1[i*n + j] = 0.0;
	result2[i*n + j] = 0.0;
	resultABC[i*n + j] = 0.0;
	resultDCB[i*n + j] = 0.0;
    }
  }
  
  printf("Multiplicando matrices de %d x %d en bloques de %d x %d\n", n, n, bs, bs);
  
  double timetick = dwalltime();
  
  MPI_Scatter(A, filas*n, MPI_DOUBLE, A, filas*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(B, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(C, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(D, filas*n, MPI_DOUBLE, D, filas*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  //Multiplicacion ABC y ABC*MAXD
  matmulblks(A, B, result1, n, bs, filas);
  matmulblks(result1, C, resultABC, n, bs, filas);
  maximoD(D, filas, n, &max);
  MPI_Allreduce(&max, &tMAX, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  numeroXMatriz(tMAX, n , filas, resultABC, P);
  
  //Multiplicacion DCB y DCB*MINA
  matmulblks(D, C, result2, n, bs, filas);
  matmulblks(result2, B, resultDCB, n, bs, filas);
  minimoA(A, filas, n, &min);
  MPI_Allreduce(&min, &tMIN, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  numeroXMatriz(tMIN, n , filas, resultDCB, P);
  
  //Calculo de R
  promedio(P, filas, n, &suma);

  MPI_Allreduce(&suma, &tSUMA, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double prom = tSUMA/(n*n);
  numeroXMatriz(prom, n, filas, P, R);
  
  MPI_Gather(R, filas*n, MPI_DOUBLE, R, filas*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  double totalTime = dwalltime() - timetick;
  
  // Validando

  int ok=0;
  for (i = 0; i<n; i++){
    for(j = 0; j<n; j++){
      if(R[i*n+j] != prom*(DN*A[i*n+j]+AN*D[i*n+j])) {
        printf("Error en %d, %d, valor: %f\n", i,j,R[i*n+j]);
      }else{ok = 1;}
    }
  }
  
  if (ok) printf("Operacion exitosa\n");

  printf("Tiempo en segundos %f\n",totalTime); 

  free(A);
  free(B);
  free(C);
  free(D);
  free(R);
  free(P);
  free(result1);
  free(result2);
  free(resultABC);
  free(resultDCB);
  
}

void ProcesoTipoB(int n, int bs, int cantProcesos){
  int i, j;
  int filas = n/cantProcesos;
  double min,max,suma,tMIN,tMAX,tSUMA;
  
  //Alocacion de porciones locales
  double *AL = (double *) malloc(sizeof(double)*filas*n);
  double *BL = (double *) malloc(sizeof(double)*n*n);
  double *CL = (double *) malloc(sizeof(double)*n*n);
  double *DL = (double *) malloc(sizeof(double)*filas*n);
  double *PL = (double *) malloc(sizeof(double)*filas*n);
  double *RL = (double *) malloc(sizeof(double)*filas*n);
  double *result1L = (double *) malloc(sizeof(double)*filas*n);
  double *result2L = (double *) malloc(sizeof(double)*filas*n);
  double *resultABCL = (double *) malloc(sizeof(double)*filas*n);
  double *resultDCBL = (double *) malloc(sizeof(double)*filas*n);
  double *V;
  
  // Inicializacion
  for (i = 0; i < filas; i++){
    for (j = 0; j < n; j++){
	PL[i*n + j] = 0;
	RL[i*n + j] = 0;
	result1L[i*n + j] = 0.0;
	result2L[i*n + j] = 0.0;
	resultABCL[i*n + j] = 0.0;
	resultDCBL[i*n + j] = 0.0;
    }
  }
  
  MPI_Scatter(V,0, MPI_DOUBLE, AL, filas*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(BL, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(CL, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(V,0, MPI_DOUBLE, DL, filas*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  

  //Multiplicacion ABC y ABC*MAXD
  matmulblks(AL, BL, result1L, n, bs, filas);
  matmulblks(result1L, CL, resultABCL, n, bs, filas);
  maximoD(DL, filas, n, &max);
  MPI_Allreduce(&max, &tMAX, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  numeroXMatriz(tMAX, n , filas, resultABCL, PL);
  
  //Multiplicacion DCB y DCB*MINA
  matmulblks(DL, CL, result2L, n, bs, filas);
  matmulblks(result2L, BL, resultDCBL, n, bs, filas);
  minimoA(AL, filas, n, &min);
  MPI_Allreduce(&min, &tMIN, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  numeroXMatriz(tMIN, n , filas, resultDCBL, PL);
  
  //Calculo de R
  promedio(PL, filas, n, &suma);
  MPI_Allreduce(&suma, &tSUMA, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  numeroXMatriz(tSUMA/(n*n), n, filas, PL, RL);
  
  MPI_Gather(RL, filas*n, MPI_DOUBLE, V, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  free(AL);
  free(BL);
  free(CL);
  free(DL);
  free(PL);
  free(RL);
  free(result1L);
  free(result2L);
  free(resultABCL);
  free(resultDCBL);
  
}

void maximoD(double *d, int filas, int n, double* MAXI){
  double max = -9999;
  for(int i = 0; i<filas; i++){
    for(int j =0; j<n; j++){
      if (d[i*n + j]>max){
        max = d[i*n + j];
      }
    }
  }
  *MAXI = max;
}

void minimoA(double *a,int filas, int n, double* MINI){
  double min = 9999;
  for(int i = 0; i<filas; i++){
    for(int j =0; j<n; j++){
      if (a[i*n + j]<min){
        min = a[i*n + j];
      }
    }
  }
  *MINI = min;
}

void numeroXMatriz(int num, int n, int filas, double *m, double * r){
  for(int i = 0; i<filas; i++){
    for(int j = 0; j<n; j++){
      r[i*n + j] += m[i*n + j] * num;
    }
  }
}

void promedio(double * P, int filas, int n, double* SUMA){
  long double suma = 0;
  for(int i = 0; i<filas; i++){
    for(int j = 0; j<n; j++){
      suma += P[i*n + j];
    }
  }
  *SUMA = suma;
}


// Multiplicación de matrices por bloques
void matmulblks(double *a, double *b, double *c, int n, int bs, int filas){
double *ablk, *bblk, *cblk;
int I, J, K;    
int i, j, k; 
 
  for(I = 0; I < filas; I += bs)
  {
    for(J = 0; J < n; J += bs)
    {
		cblk = &c[I*n + J];
		for(K = 0; K < n; K += bs)
		{
		ablk = &a[I*n + K];
		bblk = &b[J*n + K];
		
		for (i = 0; i < bs; i++)
			{
				for (j = 0; j < bs; j++)
				{
					for  (k = 0; k < bs; k++)
					{
					cblk[i*n + j] += ablk[i*n + k] * bblk[j*n + k];
					}
				}
			}
		}
    }
  }
}

int main(int argc, char *argv[]){
  MPI_Init(&argc, &argv);
  
  int id, n, bs;
  int cantidadDeProcesos;
  n = atoi(argv[1]);
  bs = atoi(argv[2]);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &cantidadDeProcesos);
  
  if (id==0) ProcesoTipoA(n, bs, cantidadDeProcesos);
  else ProcesoTipoB(n, bs, cantidadDeProcesos);
 
  MPI_Finalize();
  return 0;
}
