#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sislin.h"
#include "utils.h"
#include <likwid.h>

/* Funções auxiliares para implementar o pré-condicionador SSOR/Gauss-Seidel
   Usadas em sequência, permitem calcular z = M^−1 r = (D+ωU)^−1 * D * (D+ωL)^−1 r
*/

/* forward substitution: resolve (D + ωL) y = r (triangular inferior) */
void forward_substitution(double *D, double *L, double *r, double *y, int n, double omega) {
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        /* L is dense n x n; but non-zero entries are only where i>j */
        for (int j = 0; j < i; j++) {
            sum += omega * L[i * n + j] * y[j];
        }

        double diag = D[i * n + i];
        if (diag != 0.0)
            y[i] = (r[i] - sum) / diag;
        else {
            fprintf(stderr, "Divisão por zero na forward_substitution (D[%d,%d]==0)\n", i, i);
            exit(1);
        }
    }
}

/* backward substitution: resolve (D + ωU) z = y (triangular superior) */
void backward_substitution(double *D, double *U, double *y, double *z, int n, double omega) {
    for (int ii = n - 1; ii >= 0; ii--) {
        double sum = 0.0;
        for (int j = ii + 1; j < n; j++) {
            sum += omega * U[ii * n + j] * z[j];
        }

        double diag = D[ii * n + ii];
        if (diag != 0.0)
            z[ii] = (y[ii] - sum) / diag;
        else {
            fprintf(stderr, "Divisão por zero na backward_substitution (D[%d,%d]==0)\n", ii, ii);
            exit(1);
        }
    }
}

int main() {
    int n, k, maxit;
    double omega, epsilon;

    /* Lê parâmetros da entrada padrão e verifica erros: */
    if (scanf("%d %d %lf %d %lf", &n, &k, &omega, &maxit, &epsilon) != 5) {
        fprintf(stderr, "Erro na leitura dos parâmetros de entrada.\n");
        return 1;
    }

    if (n <= 10 || k <= 1 || k % 2 == 0) {
        fprintf(stderr, "Parâmetros inválidos: n > 10, k > 1 e k ímpar.\n");
        return 2;
    }

    /* Inicializa gerador aleatório e ponteiros */
    srandom(20252);
    likwid_markerInit();

    double *A = NULL;     /* agora A vem compactada (k*n) na v2 do sislin */
    double *b = NULL;
    double *ASP = NULL;   /* ASP é gerada pela função genSimetricaPositiva */
    double *bsp = NULL;
    double *D = NULL;
    double *L = NULL;
    double *U = NULL;
    double *M = NULL;

    double tempo_pc = 0.0;
    double tempo_dlu = 0.0;
    double tempo_precond = 0.0;

    // Gera sistema k-diagonal (A compacta) e vetor b 
    criaKDiagonal(n, k, &A, &b);

    // Gera ASP = A^T * A e bsp = A^T * b 
    genSimetricaPositiva(A, b, n, k, &ASP, &bsp, &tempo_pc);

    // Decompõe ASP em D, L, U 
    geraDLU(ASP, n, 0, &D, &L, &U, &tempo_dlu);

    /* Gera pré-condicionador M (mesma interface) */
    geraPreCond(D, L, U, omega, n, k, &M, &tempo_precond);

    /* Aloca vetores de solução */
    double *x = (double *)calloc(n, sizeof(double));
    double *x_prev = (double *)calloc(n, sizeof(double));
    if (!x || !x_prev) {
        fprintf(stderr, "Erro na alocação de memória para vetor solução.\n");
        return 3;
    }

    // Vetores auxiliares 
    double *r = (double *)calloc(n, sizeof(double));    /* resíduo */
    double *z = (double *)calloc(n, sizeof(double));    /* pré-condicionado z = M^-1 r */
    double *p = (double *)calloc(n, sizeof(double));    /* direção conjugada */
    double *Ap = (double *)calloc(n, sizeof(double));   /* ASP * p */
    if (!r || !z || !p || !Ap) {
        fprintf(stderr, "Erro na alocação de memória para vetores auxiliares.\n");
        return 4;
    }

    /* Aloca uma vez o vetor y usado em SSOR (evita malloc/free dentro do laço) */
    double *y_ssor = NULL;
    if (omega >= 1.0 && omega < 2.0) {
        y_ssor = (double *)calloc(n, sizeof(double));
        if (!y_ssor) {
            fprintf(stderr, "Erro na alocação de y_ssor.\n");
            return 5;
        }
    }

    // Inicializa r, z e p
    for (int i = 0; i < n; i++) {
        r[i] = bsp[i];

        /* Sem pré-condicionador */
        if (omega == -1) {
            z[i] = r[i];
        }
        /* Jacobi: z = D^{-1} r  (D está em D[i*n + i]) */
        else if (omega == 0.0) {
            double diag = D[i * n + i];
            if (diag == 0.0) {
                fprintf(stderr, "Divisão por zero na inicialização Jacobi D[%d,%d]==0\n", i, i);
                return 6;
            }
            z[i] = r[i] / diag;
        }
        /* SSOR / Gauss-Seidel */
        else if (omega >= 1.0 && omega < 2.0) {
            /* usa y_ssor alocado fora do loop */
            /* (D + ωL) y = r  */
            forward_substitution(D, L, r, y_ssor, n, omega);
            /* (D + ωU) z = y */
            backward_substitution(D, U, y_ssor, z, n, omega);
        }

        p[i] = z[i];
    }

    double norma = 0.0;
    double tempo_iter = 0.0;
    int iter;

    // Loop principal
    for (iter = 0; iter < maxit; iter++) {
        double temp_iter = timestamp();

        likwid_markerStartRegion("op1");
        /* Ap = ASP * p
           Observação de otimização: ASP = A^T A possui meia-banda = k-1 (porque A tem meia-banda = k/2),
           portanto ASP[i,j] == 0 para |i-j| > (k-1). Isso reduz o custo de O(n^2) para O(n*k).
        */
        int half_band_ASP = k - 1; /* meia-banda de ASP */

        for (int i = 0; i < n; i++) {
            double sum = 0.0;

            int j_start = i - half_band_ASP;
            if (j_start < 0) j_start = 0;
            int j_end = i + half_band_ASP;
            if (j_end >= n) j_end = n - 1;

            /* percorre apenas as colunas potencialmente não-nulas */
            for (int j = j_start; j <= j_end; j++) {
                sum += ASP[i * n + j] * p[j];
            }
            Ap[i] = sum;
        }

        /* alpha = (r^T z) / (p^T Ap) */
        double rz = 0.0, pAp = 0.0;
        for (int i = 0; i < n; i++) {
            rz  += r[i] * z[i];
            pAp += p[i] * Ap[i];
        }

        /* Proteção básica contra divisão por zero */
        if (pAp == 0.0) {
            fprintf(stderr, "Denominador p^T A p == 0 (possível breakdown). Iter = %d\n", iter);
            break;
        }

        double alpha = rz / pAp;

        /* x = x + alpha * p  (guarda x_prev) */
        for (int i = 0; i < n; i++) {
            x_prev[i] = x[i];
            x[i] += alpha * p[i];
        }

        /* r = r - alpha * Ap */
        for (int i = 0; i < n; i++) {
            r[i] -= alpha * Ap[i];
        }

        /* critério de parada (norma máxima entre x e x_prev) */
        norma = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = fabs(x[i] - x_prev[i]);
            if (diff > norma) norma = diff;
        }

        if (norma < epsilon)
            break;

        /* z = M^{-1} r */
        if (omega == -1) {
            for (int i = 0; i < n; i++) z[i] = r[i];
        }
        else if (omega == 0.0) {
            for (int i = 0; i < n; i++) {
                double diag = D[i * n + i];
                if (diag == 0.0) {
                    fprintf(stderr, "Divisão por zero no pré-condicionador Jacobi D[%d,%d]==0\n", i, i);
                    return 7;
                }
                z[i] = r[i] / diag;
            }
        }
        else if (omega >= 1.0 && omega < 2.0) {
            /* (D + ωL) y = r  ; (D + ωU) z = y */
            forward_substitution(D, L, r, y_ssor, n, omega);
            backward_substitution(D, U, y_ssor, z, n, omega);
        }

        /* beta = (r^T z) / (r_prev^T z_prev) */
        double rz_new = 0.0;
        for (int i = 0; i < n; i++) rz_new += r[i] * z[i];

        /* Proteção contra divisão por zero */
        double beta = 0.0;
        if (rz != 0.0)
            beta = rz_new / rz;
        else {
            fprintf(stderr, "Aviso: rz anterior == 0 (iter %d). Beta definido como 0.\n", iter);
            beta = 0.0;
        }

        /* p = z + beta * p */
        for (int i = 0; i < n; i++) {
            p[i] = z[i] + beta * p[i];
        }

        likwid_markerStopRegion("op1");
        tempo_iter += (timestamp() - temp_iter);
    }

   /* tempo médio por iteração */
    if (iter > 0) tempo_iter = tempo_iter / iter;
    else tempo_iter = 0.0;

    /* Calcula resíduo final */
    double tempo_residuo = 0.0;
    likwid_markerStartRegion("op2");
    double residuo = calcResiduoSL(A, b, x, n, k, &tempo_residuo);
    likwid_markerStopRegion("op2");

    /* Saída (mesma forma da v1) */
    printf("%d\n", n);
    for (int i = 0; i < n; i++) {
        printf("%.16g%c", x[i], (i == n - 1) ? '\n' : ' ');
    }

    printf("%.8g\n", norma);
    printf("%.16g\n", residuo);
    printf("%.8g\n", tempo_precond);
    printf("%.8g\n", tempo_iter);
    printf("%.8g\n", tempo_residuo);

    /* Libera memória */
    free(A); free(b); free(ASP); free(D); free(L); free(U); free(M);
    free(x); free(x_prev); free(r); free(z); free(p); free(Ap); free(bsp);
    if (y_ssor) 
        free(y_ssor);

    likwid_markerClose();
    return 0;
}
