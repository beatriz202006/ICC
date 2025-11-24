#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"


// Funções auxiliares para acessar matriz k-diagonal compacta

static inline int idxDiag(int i, int j, int k) {
    return j - i + k/2;
}

static inline double getA(double *A, int n, int k, int i, int j) {
    int d = idxDiag(i, j, k);
    if (d < 0 || d >= k) 
      return 0.0;

    return A[d*n + i];
}

static inline void setA(double *A, int n, int k, int i, int j, double val) {
    int d = idxDiag(i, j, k);
    if (d < 0 || d >= k) 
      return;
    
    A[d*n + i] = val;
}

// Funções auxiliares originais

static inline double generateRandomA(unsigned int i, unsigned int j, unsigned int k) {
    static double invRandMax = 1.0 / (double)RAND_MAX;
    return ((i==j) ? (double)(k<<1) : 1.0) * (double)random() * invRandMax;
}

static inline double generateRandomB(unsigned int k) {
    static double invRandMax = 1.0 / (double)RAND_MAX;
    return (double)(k<<2) * (double)random() * invRandMax;
}

// Criação da matriz k-diagonal compacta

void criaKDiagonal(int n, int k, double **A, double **B)
{
    // Agora A tem k*n posições (apenas diagonais)
    *A = (double *)calloc(k * n, sizeof(double));
    *B = (double *)calloc(n, sizeof(double));

    if (!(*A) || !(*B)) {
        fprintf(stderr, "Erro na alocação de memória para matriz A ou vetor B.\n");
        exit(1);
    }

    // Preenche somente diagonais não nulas
    for (int i = 0; i < n; i++) {
        for (int j = i - k/2; j <= i + k/2; j++) {
            if (j >= 0 && j < n) {
                double val = generateRandomA(i, j, k);
                setA(*A, n, k, i, j, val);
            }
        }
    }

    for (int i = 0; i < n; i++)
        (*B)[i] = generateRandomB(k);
}


// ASP = A^T * A (ambas k-diagonais)


void genSimetricaPositiva(double *A, double *b, int n, int k, double **ASP, double **bsp, double *tempo)
{
    *tempo = timestamp();

    *ASP = (double *)calloc(n * n, sizeof(double));
    *bsp = (double *)calloc(n, sizeof(double));

    if (!(*ASP) || !(*bsp)) {
        fprintf(stderr, "Erro na alocação de memória para matriz ASP ou vetor BSP.\n");
        exit(1);
    }

    // ASP = A^T * A
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int l = 0; l < n; l++) {
                sum += getA(A, n, k, l, i) * getA(A, n, k, l, j);
            }
            (*ASP)[i*n + j] = sum;
        }
    }

    // bsp = A^T * b
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int l = 0; l < n; l++) {
            sum += getA(A, n, k, l, i) * b[l];
        }
        (*bsp)[i] = sum;
    }

    *tempo = timestamp() - *tempo;
}


// Decomposição DLU apenas reorganiza valores existentes
void geraDLU(double *A, int n, int k_unused, double **D, double **L, double **U, double *tempo)
{
    *tempo = timestamp();

    *D = calloc(n*n, sizeof(double));
    *L = calloc(n*n, sizeof(double));
    *U = calloc(n*n, sizeof(double));

    if (!(*D) || !(*L) || !(*U)) {
        fprintf(stderr, "Erro na alocação de memória para D, L e U.\n");
        exit(1);
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {

            double val = A[i*n + j];   

            if (i == j) (*D)[i*n + j] = val;      // D
            else if (i > j) (*L)[i*n + j] = val;  // L
            else (*U)[i*n + j] = val;             // U
        }
    }

    *tempo = timestamp() - *tempo;
}

// Pré-condicionador igual ao v1
void geraPreCond(double *D, double *L, double *U, double w, int n, int k,
                 double **M, double *tempo)
{
    *tempo = timestamp();

    *M = (double *)calloc(n * n, sizeof(double));

    if (!(*M)) {
        fprintf(stderr, "Erro na alocação de memória para M.\n");
        exit(1);
    }

    if (w == -1) {
        for (int i = 0; i < n; i++)
            (*M)[i*n + i] = 1.0;
    }
    else if (w == 0.0) {
        for (int i = 0; i < n; i++)
            (*M)[i*n + i] = D[i*n + i];
    }
    else if (w >= 1.0 && w < 2.0) {
        for (int i = 0; i < n; i++)
            (*M)[i*n + i] = 1.0;
    }
    else {
        fprintf(stderr, "Valor de w inválido.\n");
        exit(1);
    }

    *tempo = timestamp() - *tempo;
}

// Cálculo do resíduo usando matriz k-diagonal compacta
double calcResiduoSL(double *A, double *b, double *X, int n, int k, double *tempo)
{
    *tempo = timestamp();

    double *r = calloc(n, sizeof(double));
    if (!r) {
        fprintf(stderr, "Erro na alocação de r.\n");
        exit(1);
    }

    for (int i = 0; i < n; i++) {
        double Ax_i = 0.0;

        for (int j = i - k/2; j <= i + k/2; j++) {
            if (j >= 0 && j < n) {
                Ax_i += getA(A, n, k, i, j) * X[j];
            }
        }

        r[i] = b[i] - Ax_i;
    }

    double norma = 0.0;
    for (int i = 0; i < n; i++)
        norma += r[i] * r[i];

    free(r);
    *tempo = timestamp() - *tempo;

    return sqrt(norma);
}
