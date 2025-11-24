#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sislin.h"
#include "utils.h"
#include <likwid.h>

/* Funções auxiliares para implementar o pré-condicionador SSOR/Gauss-Seidel
Usadas em sequência, permitem calcular z =M^−1r = (D+ωU)^−1*D(D+ωL)^−1r */

/* Função para forward substitution: resolve (D + ωL)y = r (sistema triangular inferior)
y é o vetor de saída que recebe a solução do sistema triangular inferior */
void forward_substitution(double *D, double *L, double *r, double *y, int n, double omega) {
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += omega * L[i * n + j] * y[j]; //acumula a soma de omega * L[i][j] * y[j]
        }

        if ((D[i * n + i]) != 0)
            y[i] = (r[i] - sum) / (D[i * n + i]);
        
        else {
            printf("Divisão por zero!");
            exit(1);
        }
    }
} 

/* Função para backward substitution: resolve (D + ωU)z = y (sistema triangular superior)
z é o vetor de saída que recebe a solução do sistema triangular superior */
void backward_substitution(double *D, double *U, double *y, double *z, int n, double omega) {
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += omega * U[i * n + j] * z[j]; //acumula a soma de omega * U[i][j] * z[j]
        }

        if ((D[i * n + i]) != 0)
            z[i] = (y[i] - sum) / (D[i * n + i]);
        
        else {
            printf("Divisão por zero!");
            exit(1);
        }
    }
}

int main() {
    int n, k, maxit;
    double omega, epsilon;

    //Lê parâmetros da entrada padrão e verifica erros:
    if (scanf("%d %d %lf %d %lf", &n, &k, &omega, &maxit, &epsilon) != 5) {
        fprintf(stderr, "Erro na leitura dos parâmetros de entrada.\n");
        return 1;
    }

    if (n <= 10 || k <= 1 || k % 2 == 0) {
        fprintf(stderr, "Parâmetros inválidos: n > 10, k > 1 e k ímpar.\n");
        return 2;
    }

    //Inicializa o gerador de números aleatórios, ponteiros de matrizes e vetores:
    srandom(20252);
    likwid_markerInit();

    double *A = NULL;
    double *b = NULL;
    double *ASP = NULL;
    double *bsp = NULL;
    double *D = NULL;
    double *L = NULL;
    double *U = NULL;
    double *M = NULL;

    double tempo_pc = 0.0;
    double tempo_dlu = 0.0;
    double tempo_precond = 0.0;

    //Gera sistema k-diagonal:
    criaKDiagonal(n, k, &A, &b);

    //Gera matriz simétrica positiva e vetor bsp:
    genSimetricaPositiva(A, b, n, k, &ASP, &bsp, &tempo_pc);

    //Decompõe ASP em D, L e U:
    geraDLU(ASP, n, k, &D, &L, &U, &tempo_dlu);

    //Gera pré-condicionador M:
    geraPreCond(D, L, U, omega, n, k, &M, &tempo_precond);

    //Aloca vetor x de solução e o vetor x_prev que mantém valores da iteração anterior
    double *x = (double *)calloc(n, sizeof(double));
    double *x_prev = (double *)calloc(n, sizeof(double));

    //Verifica se houve erro na alocação:
    if (!x || !x_prev) {
        fprintf(stderr, "Erro na alocação de memória para vetor solução.\n");
        return 3;
    }

    //Vetores auxiliares:
    double *r = (double *)calloc(n, sizeof(double));    //vetor do resíduo
    double *z = (double *)calloc(n, sizeof(double));    //vetor pré-condicionado: z = M^-1 r
    double *p = (double *)calloc(n, sizeof(double));    //vetor direção do gradiente conjugado
    double *Ap = (double *)calloc(n, sizeof(double));   // A*p p/ calculo do alpha

    //Verifica se houve erro na alocação:
    if (!r || !z || !p || !Ap) {
        fprintf(stderr, "Erro na alocação de memória para vetores auxiliares.\n");
        return 4;
    }

    //Inicializa r, z e p:
    for (int i = 0; i < n; i++) {
        r[i] = bsp[i];
        
        // Sem pré-condicionador
        if (omega == -1) {
            z[i] = r[i]; 
        }

        // Com pré-condicionador de Jacobi
        else if (omega == 0.0) {
            z[i] = r[i] / D[i * n + i]; 
        }

        //SSOR/ Gauss-Seidel
        else if (omega >= 1.0 && omega < 2.0) {

            //Aloca vetor y para forward substitution:
            double *y = (double *)calloc(n, sizeof(double));
            forward_substitution(D, L, r, y, n, omega);
            backward_substitution(D, U, y, z, n, omega);
            free(y);
        }

        p[i] = z[i];
    }

    double norma = 0.0;
    double tempo_iter = 0.0;
    int iter;

    //Loop principal (executa até maxit ou até convergência)
    for (iter = 0; iter < maxit; iter++) {
        double temp_iter = timestamp();

        likwid_markerStartRegion("op1");

        //Ap = ASP * p
        for (int i = 0; i < n; i++) {
            Ap[i] = 0.0;
            for (int j = 0; j < n; j++) {
                Ap[i] += ASP[i * n + j] * p[j]; //calcula Ap = ASP[i][j] * p[j] => produto matriz-vetor
            }
        }

        // alpha = (r^T z) / (p^T Ap) => formula
        double rz = 0.0, pAp = 0.0;
        for (int i = 0; i < n; i++) {
            rz += r[i] * z[i];    //numerador de alpha
            pAp += p[i] * Ap[i];  //denominador de alpha
        }

        double alpha = rz / pAp;

        //x = x + alpha * p => multiplica a direção p pelo tamanho do passo α e atualiza a solução x
        for (int i = 0; i < n; i++) {
            x_prev[i] = x[i];
            x[i] += alpha * p[i];
        }

        // Atualiza resíduo: r = r - alpha * Ap
        for (int i = 0; i < n; i++) {
            r[i] -= alpha * Ap[i];
        }

        //Verifica critério de parada (norma máxima)
        norma = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = fabs(x[i] - x_prev[i]);
            if (diff > norma)
                norma = diff;
        }

        // Verifica se convergiu
        if (norma < epsilon)
            break;
        
        //z = M^-1 r (Jacobi ou identidade)
        for (int i = 0; i < n; i++) {

            //Sem pré-condicionador
            if (omega == -1) {
                z[i] = r[i];
            }

            //Com pré-condicionador de Jacobi
            else if (omega == 0.0) {
                z[i] = r[i] / D[i * n + i];
            }

            //SSOR/ Gauss-Seidel
            else if (omega >= 1.0 && omega < 2.0) {
                double *y = (double *)calloc(n, sizeof(double));
                forward_substitution(D, L, r, y, n, omega);
                backward_substitution(D, U, y, z, n, omega);
                free(y);
            }
        }

        //beta = (r^T * z) / (r_prev^T * z_prev) => beta ajusta a contribuição da direção antiga p_k na nova direção
        double rz_new = 0.0;
        for (int i = 0; i < n; i++) {
            rz_new += r[i] * z[i];
        }

        double beta = rz_new / rz;

        //p = z + beta *  p => atualiza a nova direção de busca p
        for (int i = 0; i < n; i++) {
            p[i] = z[i] + beta * p[i];
        }

        likwid_markerStopRegion("op1");
        tempo_iter += (timestamp() - temp_iter);
    }

    //Calcula tempo médio por iteração:
    if (iter > 0)
        tempo_iter = tempo_iter / iter;
    else
        tempo_iter = 0.0;
    
    //Calcula resíduo final:
    double tempo_residuo = 0.0;
    likwid_markerStartRegion("op2");
    double residuo = calcResiduoSL(A, b, x, n, k, &tempo_residuo);
    likwid_markerStopRegion("op2");

    //Imprime resultados:
    printf("%d\n", n);
    for (int i = 0; i < n; i++) {
        printf("%.16g%c", x[i], (i == n - 1) ? '\n' : ' ');
    }

    printf("%.8g\n", norma);
    printf("%.16g\n", residuo);
    printf("%.8g\n", tempo_precond);
    printf("%.8g\n", tempo_iter);
    printf("%.8g\n", tempo_residuo);

    // Libera memória
    free(A); free(b); free(ASP); free(D); free(L); free(U); free(M);
    free(x); free(x_prev); free(r); free(z); free(p); free(Ap); free(bsp);

    likwid_markerClose();
    return 0;
}