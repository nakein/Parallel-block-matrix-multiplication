#include<stdio.h>
#include<stdlib.h>  
#include <sys/time.h>
#include <omp.h>
#include <time.h>

double dwalltime(){
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

double maximoD(double *d, int n){
  int maxi = -9999;
  #pragma omp parallel for reduction(max:maxi)
  for(int i = 0; i<n; i++){
    for(int j =0; j<n; j++){
      if (d[i*n + j]>maxi){
        maxi = d[i*n + j];
      }
    }
  }
  return maxi;
}

double minimoA(double *a, int n){
  int mini = 9999;
  #pragma omp parallel for reduction(min:mini) 
  for(int i = 0; i<n; i++){
    for(int j =0; j<n; j++){
      if (a[i*n + j]<mini){
        mini = a[i*n + j];
      }
    }
  }
  return mini;
}

void numeroXMatriz(int num, int n, double *m, double * r){
  #pragma omp parallel for
  for(int i = 0; i<n; i++){
    for(int j = 0; j<n; j++){
      r[i*n + j] += m[i*n + j] * num;
    }
  }
}

double promedio(double * P, int n){
  long double suma = 0;
  #pragma omp parallel for reduction(+:suma)
  for(int i = 0; i<n; i++){
    for(int j = 0; j<n; j++){
      suma += P[i*n + j];
    }
  }
  
  return suma/(n*n);
}


// Multiplicación de matrices por bloques
void matmulblks(double *a, double *b, double *c, int n, int bs){
double *ablk, *bblk, *cblk;
int I, J, K;    
int i, j, k; 
 
  #pragma omp parallel for private(I,J,K, ablk, bblk, cblk, i,j,k)
  for(I = 0; I < n; I += bs)
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
  double *A, *B, *C, *D, *R, *P, *result1, *result2, *resultABC, *resultDCB;
  int n, bs, i, j, T;

  double timetick;
  
  T = atoi(argv[3]);	
  
    // Chequeo de parámetros
	if ( (argc != 4) || ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((n % bs) != 0)){
		printf("Error en los parámetros. Usar: ./%s N BS (N debe ser multiplo de BS)\n", argv[0]);
		exit(1);
	}
	
  //Set numero de hilos
  omp_set_num_threads(T);

  // Alocar  
	A = (double *) malloc(n*n*sizeof(double));
	B = (double *) malloc(n*n*sizeof(double));
	C = (double *) malloc(n*n*sizeof(double));
	D = (double *) malloc(n*n*sizeof(double));
	R = (double *) malloc(n*n*sizeof(double));
	P = (double *) malloc(n*n*sizeof(double));
	result1 = (double *) malloc(n*n*sizeof(double));
	result2 = (double *) malloc(n*n*sizeof(double));
	resultABC = (double *) malloc(n*n*sizeof(double));
	resultDCB = (double *) malloc(n*n*sizeof(double));

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
			R[i*n + j] = 0.0;
			P[i*n + j] = 0.0;
			result1[i*n + j] = 0.0;
			result2[i*n + j] = 0.0;
			resultABC[i*n + j] = 0.0;
			resultDCB[i*n + j] = 0.0;
		}
	}


	printf("Multiplicando matrices de %d x %d en bloques de %d x %d\n", n, n, bs, bs);
  
  timetick = dwalltime();  
    
  //Multiplicacion ABC y ABC*MAXD
  matmulblks(A, B, result1, n, bs);
  matmulblks(result1, C, resultABC, n, bs);
  double maxD = maximoD(D, n);
  numeroXMatriz(maxD, n , resultABC, P);
  
  //Multiplicacion DCB y DCB*MINA
  matmulblks(D, C, result2, n, bs);
  matmulblks(result2, B, resultDCB, n, bs);
  double minA = minimoA(A, n);
  numeroXMatriz(minA, n , resultDCB, P);
  
  //Calculo de R
  double prom = promedio(P,n);
  numeroXMatriz(prom, n, P, R);
  
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
 
	return 0;
}
