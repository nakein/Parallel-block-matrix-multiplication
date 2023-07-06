#include<stdio.h>
#include<stdlib.h>  
#include <sys/time.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

// Prototipos
void matmulblks(double *, double *, double *, int, int, int, int);
void maximoD(double *, double *, int, int, int);
void numeroXMatriz(double, int, double *, double * , int, int);
void minimoA(double *, double *, int, int, int);
void promedio(double *, double *, int, int, int T);

// Compartidas
pthread_barrier_t* barrera;
double *A, *B, *C, *D, *R, *P, *result1, *result2, *resultABC, *resultDCB;
double *max, *min, *prom;
int n, bs, T;

double dwalltime(){
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

void * funcion(void *arg){
  int tid = *(int*)arg;
  
  //Primera parte de P
  matmulblks(A, B, result1, n, bs, tid, T);
  matmulblks(result1, C, resultABC, n, bs, tid, T);
  maximoD(max, D, n, tid, T);
  numeroXMatriz(max[0], n , resultABC, P, tid, T);
  
  //Segunda parte de P
  matmulblks(D, C, result2, n, bs, tid, T);
  matmulblks(result2, B, resultDCB, n, bs, tid, T);
  minimoA(min, A, n, tid, T);
  numeroXMatriz(min[0], n , resultDCB, P, tid, T);
  
  //Resultado R
  promedio(prom, P, n, tid, T);
  numeroXMatriz(prom[0], n, P, R, tid, T);
  
  pthread_exit(NULL);
}

void maximoD(double *max, double *d, int n, int tid, int T){
  max[tid] = -9999;
  int porciones, porciones2;
  int limitebajo = (n*tid)/T;
  int limitealto = (n*(tid+1))/T;
  for(int i = limitebajo; i<limitealto; i++){
    for(int j =0; j<n; j++){
      if (d[i*n + j]>max[tid]){
        max[tid] = d[i*n + j];
      }
    }
  }
  pthread_barrier_wait(&barrera[0]);  
  
  for (int i=1, j=0; i<=log2(T); i++, j++){
        porciones = pow(2,i);
	porciones2 = pow(2,j);
        if (tid % porciones != 0) break;
        
        if(max[tid] < max[tid+porciones2]){
          max[tid] = max[tid+porciones2];
          max[tid+porciones2] = 0;
        }

        pthread_barrier_wait(&barrera[i]);
  }
  pthread_barrier_wait(&barrera[0]); 
}

void minimoA(double *min, double *a, int n, int tid, int T){
  min[tid] = 9999;
  int limitebajo = (n*tid)/T;
  int limitealto = (n*(tid+1))/T;
  int porciones, porciones2;
  for(int i = limitebajo; i<limitealto; i++){
    for(int j =0; j<n; j++){
      if (a[i*n + j]<min[tid]){
        min[tid] = a[i*n + j];
      }
    }
  }

  pthread_barrier_wait(&barrera[0]);  
  
  
  for (int i=1, j=0; i<=log2(T); i++, j++){
        porciones = pow(2,i);
	porciones2 = pow(2,j);
        if (tid % porciones != 0) break;
        
        if(min[tid] > min[tid+porciones2]){
          min[tid] = min[tid+porciones2];
          min[tid+porciones2] = 0;
        }

        pthread_barrier_wait(&barrera[i]);
  }
  pthread_barrier_wait(&barrera[0]); 
}

void numeroXMatriz(double num, int n, double *m, double * r, int tid, int T){
  int limitebajo = (n*tid)/T;
  int limitealto = (n*(tid+1))/T;
  for(int i = limitebajo; i<limitealto; i++){
    for(int j = 0; j<n; j++){
      r[i*n + j] += m[i*n + j] * num;
    }
  }
}

void promedio(double *prom, double * P, int n, int tid, int T){
  prom[tid]=0;
  int limitebajo = (n*tid)/T;
  int limitealto = (n*(tid+1))/T;
  int porciones, porciones2;
  for(int i = limitebajo; i<limitealto; i++){
    for(int j = 0; j<n; j++){
      prom[tid] += P[i*n + j];
    }
  }
  
  pthread_barrier_wait(&barrera[0]); 
  
  
  for (int i=1, j=0; i<=log2(T); i++, j++){
        porciones = pow(2,i);
	porciones2 = pow(2,j);
        if (tid % porciones != 0) break;
         
         prom[tid] = prom[tid] + prom[tid+porciones2];
         prom[tid+porciones2] = 0;

        pthread_barrier_wait(&barrera[i]);
  }
  
  if(tid == 0){
  
  prom[tid] = prom[tid]/(n*n);}
  
  pthread_barrier_wait(&barrera[0]); 

}


// Multiplicación de matrices por bloques
void matmulblks(double *a, double *b, double *c, int n, int bs, int tid, int T){
double *ablk, *bblk, *cblk;
int I, J, K;    
int i, j, k; 
int limitebajo = (n*tid)/T;
int limitealto = (n*(tid+1))/T;

  for(I = limitebajo; I < limitealto; I += bs)
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
  int i, j;
  T = atoi(argv[3]);
  pthread_t misThreads[T];
  int threads_ids[T];
  
  barrera = (pthread_barrier_t*) malloc(T * sizeof(pthread_barrier_t));
  for(i=0; i<=log2(T); i++) pthread_barrier_init(&barrera[i], NULL, T/pow(2,i));

  double timetick;

  // Chequeo de parámetros
	if ( (argc != 4) || ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((n % bs) != 0)){
		printf("Error en los parámetros. Usar: ./%s N BS (N debe ser multiplo de BS)\n", argv[0]);
		exit(1);
	}

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

	min = malloc(T*sizeof(double));
	max = malloc(T*sizeof(double));
	prom = malloc(T*sizeof(double));

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
  
  for(int id=0; id<T; id++){
    threads_ids[id] = id;
    pthread_create(&misThreads[id], NULL, &funcion,(void*)&threads_ids[id]);
  }
  
  for(int id=0; id<T; id++){
    pthread_join(misThreads[id], NULL);
  }
  
  double totalTime = dwalltime() - timetick;

  // Validando
  int ok=0;
  for (i = 0; i<n; i++){
    for(j = 0; j<n; j++){
      if(R[i*n+j] != prom[0]*(DN*A[i*n+j]+AN*D[i*n+j])) {
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
	
	for(i=0; i<=log2(T); i++) pthread_barrier_destroy(&barrera[i]);
	free(barrera);

	return 0;
}
