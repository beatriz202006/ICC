#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main() {

  int n,k,it;
  double **M, *Y;
  double t, norma;

  while ( (k=scanf("%d\n", &n)) != EOF) {
    if (k != 1) {
      printf("Formato de saída INCORRETO [n]\n");
      exit(1);
    }
  
    M = (double **) malloc(sizeof(double *) * n);
    for (int i = 0; i < n; ++i)
      M[i] = (double *) malloc(sizeof(double) * n+1);
    
    for (int i =0; i < n; ++i)
      for (int j = 0; j < n+1; ++j)
	if (scanf("%lf", &(M[i][j])) != 1) {
	  printf("Formato de saída INCORRETO [A | B]\n");
	  exit(1);
	}

    Y = (double *) malloc(sizeof(double) * n);
    for (int j = 0; j < n; ++j)
      if (scanf("%lf", &(Y[j])) != 1) {
	printf("Formato de saída INCORRETO [Y]\n");
	exit(1);
      }
  
    if (scanf("%d\n", &it) != 1) {
      printf("Formato de saída INCORRETO [iterações]\n");
      exit(1);
    }

    if (scanf("%lf\n", &norma) != 1) {
      printf("Formato de saída INCORRETO [normaL2]\n");
      exit(1);
    }

    if (scanf("%lf\n", &t) != 1) {
      printf("Formato de saída INCORRETO [tempo]\n");
      exit(1);
    }

    printf("n = %d\n", n);
    printf("A | B:\n");
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n+1; ++j)
	printf("%23.15e ", M[i][j]);
      printf("\n");
    }

    printf("Y:\n");
    for (int j = 0; j < n; ++j)
      printf("%23.15e ", Y[j]);
    printf("\n");

    printf("iterações = %d\n", it);
    printf("norma L2 = %23.15e\n", norma);
    printf("t = %16.8e\n", t);
    printf("Formato de saída correto!\n");
  
    for (int i = 0; i < n; ++i)
      free(M[i]);
    free(M);
    free(Y);
  }
}
