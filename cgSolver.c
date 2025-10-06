#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sislin.h"
#include "utils.h"

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

    //Inicializa o gerador de números aleatórios:
    srandom(20252);

    //Gera sistema k-diagonal:
    double *A = NULL;
    double *b = NULL;
    criaKDiagonal(n, k, &A, &b);

    //Gera matriz simétrica positiva e vetor bsp:
    double *ASP = NULL;
    double *bsp = NULL;
    double tempo_pc = 0.0;
    genSimetricaPositiva(A, b, n, k, &ASP, &bsp, &tempo_pc);

    //Decompõe ASP em D, L e U:
    double *D = NULL;
    double *L = NULL;
    double *U = NULL;
    double tempo_dlu = 0.0;
    geraDLU(ASP, n, k, &D, &L, &U, &tempo_dlu);

    //Gera pré-condicionador M:
    double *M = NULL;
    double tempo_precond = 0.0;
    geraPreCond(D, L, U, omega, n, k, &M, &tempo_precond);

    //Inicializa a solução X com zeros:
    double *x = (double *)calloc(n, sizeof(double));
    double *x_prev = (double *)calloc(n, sizeof(double));

    //Verifica se houve erro na alocação:
    if (!x || !x_prev) {
        fprintf(stderr, "Erro na alocação de memória para vetor solução.\n");
        return 3;
    }

    //Vetores auxiliares:
    double *r = (double *)calloc(n, sizeof(double));
    double *z = (double *)calloc(n, sizeof(double));
    double *p = (double *)calloc(n, sizeof(double));
    double *Ap = (double *)calloc(n, sizeof(double));

    //Verifica se houve erro na alocação:
    if (!r || !z || !p || !Ap) {
        fprintf(stderr, "Erro na alocação de memória para vetores auxiliares.\n");
        return 4;
    }

    //Inicializa r = b - ASP * x ( x= 0 => r = b)
    for (int i = 0; i < n; i++) {
        r[i] = b[i];
        
        // Sem pré-condicionador
        if (omega == -1) {
            z[i] = r[i]; 
        }

        // Com pré-condicionador de Jacobi
        else if (omega == 0.0) {
            z[i] = r[i] / D[i * n + i]; 
        }

        p[i] = z[i];
    }

    double norma = 0.0;
    double tempo_iter = 0.0;
    int iter;

    for (iter = 0; iter < maxit; iter++) {
        double temp_iter = timestamp();

        //Ap = ASP * p
        for (int i = 0; i < n; i++) {
            Ap[i] = 0.0;
            for (int j = 0; j < n; j++) {
                Ap[i] += ASP[i * n + j] * p[j];
            }
        }

        // alpha = (r^T z) / (p^T Ap)
        double rz = 0.0, pAp = 0.0;
        for (int i = 0; i < n; i++) {
            rz += r[i] * z[i];
            pAp += p[i] * Ap[i];
        }

        double alpha = rz / pAp;

        //x = x + alpha * p
        for (int i = 0; i < n; i++) {
            x_prev[i] = x[i];
            x[i] += alpha * p[i];
        }

        // r = r - alpha * Ap
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
        }

        //beta = r^T z) / (r_prev^T z_prev)
        double rz_new = 0.0, rz_old = 0.0;
        for (int i = 0; i < n; i++) {
            rz_new += r[i] * z[i];
            rz_old += r[i] * z[i];
        }

        double beta = rz_new / rz;

        //p = z + beta *  p
        for (int i = 0; i < n; i++) {
            p[i] = z[i] + beta * p[i];
        }

        tempo_iter += (timestamp() - tempo_iter);
    }

    //Calcula tempo médio por iteração:
    if (iter > 0)
        tempo_iter = tempo_iter / iter;
    else
        tempo_iter = 0.0;
    
    //Calcula resíduo final:
    double tempo_residuo = 0.0;
    double residuo = calcResiduoSL(ASP, bsp, x, n, k, &tempo_residuo);

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

    return 0;
}