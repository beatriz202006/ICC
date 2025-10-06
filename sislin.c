#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k );
static inline real_t generateRandomB( unsigned int k );

/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return ( (i==j) ? (real_t)(k<<1) : 1.0 )  * (real_t)random() * invRandMax;
}

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomB( unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return (real_t)(k<<2) * (real_t)random() * invRandMax;
}


/* Cria matriz 'A' k-diagonal e Termos independentes B */
void criaKDiagonal(int n, int k, real_t **A, real_t **B)
{
  // Aloca matriz A (n X n) e vetor B (n)
  *A = (real_t *)calloc(n * n, sizeof(real_t));
  *B = (real_t *)calloc(n, sizeof(real_t));

  //Verifica se houve erro na alocação:
  if (!(*A) || !(*B)) {
    fprintf(stderr, "Erro na alocação de memória para matriz A ou vetor B.\n");
    exit(1);
  }

  //Preenche a matriz A ( k-diagonal):
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      //Preenche se estiver em uma das k diagonais:
      if (abs(i-j) <= k/2) {
        (*A)[i * n + j] = generateRandomA(i, j, k);
      }
    }
  }

  //Preenche o vetor B:
  for (int i = 0; i < n; i++) {
    (*B)[i] = generateRandomB(k);
  }
}

/* Gera matriz simetrica positiva */
void genSimetricaPositiva(double *A, double *b, int n, int k, 
			  double **ASP, double **bsp, double *tempo)
{
  *tempo = timestamp();

  //Aloca matriz ASP (n x n) :
  *ASP = (double *)calloc(n * n, sizeof(double));

  //Aloca vetor bsp (n):
  *bsp = (double *)calloc(n, sizeof(double));

  //Verifica se houve erro na alocação:
  if (!(*ASP) || !(*bsp)) {
    fprintf(stderr, "Erro na alocação de memória para matriz ASP ou vetor BSP.\n");
    exit(1);
  }

  //Calcula ASP = A * A ^T:
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0.0;
      for (int l = 0; l < n; l++) {
        sum += A[i * n + l] * A[j * n + l];
      }     

      (*ASP)[i * n + j] = sum;
    }
  }

  //Copia b para bsp:
  for (int i = 0; i < n; i++) {
    (*bsp)[i] = b[i];
  }

  *tempo = timestamp() - *tempo;
}
  
/* Decompõe A em D, L e U */
void geraDLU (double *A, int n, int k,
	      double **D, double **L, double **U, double *tempo)
{
  *tempo = timestamp();
  //Aloca matrizes D, L e U (n x n):
  *D = (double *)calloc(n * n, sizeof(double));
  *L = (double *)calloc(n * n, sizeof(double));
  *U = (double *)calloc(n * n, sizeof(double));

  //Verifica se houve erro na alocação:
  if (!(*D) || !(*L) || !(*U)) {
    fprintf(stderr, "Erro na alocação de memória para matrizes D, L ou U.\n");
    exit(1);
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      
      double val = A[i * n + j];
      //Se estiver na diagonal principal:
      if (i == j) { 
        (*D)[i * n + j] = val;
      }

      //Se estiver na parte inferior:
      else if (i > j) {
        (*L)[i * n + j] = val;
      }

      //Se estiver na parte superior:
      else if ( i < j) {
        (*U)[i * n + j] = val;
      }
    }
  }

  *tempo = timestamp() - *tempo;
}
  
/**
 * Devolve matriz M⁻¹
 *
 */
void geraPreCond(double *D, double *L, double *U, double w, int n, int k,
		 double **M, double *tempo)
{
  *tempo = timestamp();
  //Aloca matriz M (n * n):
  *M = (double *)calloc(n * n, sizeof(double));

  //Verifica se houve erro na alocação:
  if ((!*M)) {
    fprintf(stderr, "Erro na alocação de memória para matriz M.\n");
    exit(1);
  }

  //Sem pré-condicionador: M = I (identidade)
  if (w == -1) {
    for (int i = 0; i < n; i++) {
      (*M)[i * n + i] = 1.0;
    }
  }

  //Com pré-condicionador de Jacobi: M = D
  else if (w == 0.0) {
    for (int i = 0; i < n; i++) {
      (*M)[i * n + i] = D[i * n + i];
    }
  }

  //SSOR/ Gauss-Seidel: M = (D + ωL)D^-1 (D + ωU )
  else if (w >= 1.0 && w < 2.0){
    for (int i = 0; i < n; i++) {
    (*M)[i * n + i] = 1.0;
    }
  }

  else {
    fprintf(stderr, "Valor de w inválido para pré-condicionador.\n");
    exit(1);
  }

  *tempo = timestamp() - *tempo;
}

double calcResiduoSL (double *A, double *b, double *X,
		      int n, int k, double *tempo)
{
  *tempo = timestamp();

  //Aloca vetor r (n):
  double *r = calloc(n, sizeof(double));

  //Verifica se houve erro na alocação:
  if (!r) {
    fprintf(stderr, "Erro na alocação de memória para vetor resíduo.\n");
    exit(1);
  }

  //Calcula r = b - A*x:
  for (int i = 0; i < n; i++) {
    double Ax_i = 0.0;
    for (int j = 0; j < n; j++) {
      Ax_i += A[i * n + j] * X[j];
    }

    r[i] = b[i] - Ax_i;
  }

  //Calcula norma L2 do resíduo:
  double norma = 0.0;
  for (int i = 0; i < n; i++) {
    norma += r[i] * r[i];
  }

  norma  = sqrt(norma);
  free(r);
  *tempo = timestamp() - *tempo;

  return norma;
}


