/*Beatriz Pontes Camargo 
GRR 20242966
*/
#include <stdio.h>
#include <stdlib.h>
#include <fenv.h>
#include "edo.h"
#include "gauss_seidel.h"
#include "utils.h"
#include <likwid.h>


int main() {
    // Força arredondamento para baixo
    fesetround(FE_DOWNWARD);

    //Inicializa os marcadores do LIKWID
    LIKWID_MARKER_INIT;

    //Declaração das variáveis para os parâmetros do problema
    int n;
    real_t a, b, ya, yb, p, q;
    EDo edo;
    real_t r1, r2, r3, r4;

    // Lê dados iniciais                        
    if (scanf("%d", &n) != 1) return 1;             // Número de pontos
    if (scanf("%lf %lf", &a, &b) != 2) return 1;    // Intervalo [a, b]
    if (scanf("%lf %lf", &ya, &yb) != 2) return 1;  // Valores de contorno ya, yb
    if (scanf("%lf %lf", &p, &q) != 2) return 1;    // Parâmetros p, q

    // Preenche a estrutura EDO com os dados lidos
    edo.n = n;
    edo.a = a;
    edo.b = b;
    edo.ya = ya;
    edo.yb = yb;
    edo.p = p;
    edo.q = q;
    
    // Para cada linha de coeficientes r1 r2 r3 r4
    // Cada conjunto representa uma nova EDO a ser resolvida
    while (scanf("%lf %lf %lf %lf", &r1, &r2, &r3, &r4) == 4) {

        //Atualiza os coeficientes na estrutura EDO
        edo.r1 = r1;
        edo.r2 = r2;
        edo.r3 = r3;
        edo.r4 = r4;

        // Gera sistema tridiagonal correspondente
        Tridiag *sl = genTridiag(&edo);

        // Imprime o sistema gerado (matriz e termos independentes)
        prnEDOsl(&edo);

        //Aloca vetor solução
        real_t *Y = (real_t *)calloc(n, sizeof(real_t));
        int maxIt = 100;
        real_t tol = 1e-5;
        real_t normaL2;
        int numIt;

        // Marca o início da região a ser medida pelo LIKWID
        LIKWID_MARKER_START("GaussSeidel");
        rtime_t t0 = timestamp(); // marca o tempo inicial

        // Resolve o sistema usando Gauss-Seidel
        gaussSeidel(sl, Y, maxIt, tol, &normaL2, &numIt); 
        rtime_t t1 = timestamp() - t0; // Tempo gasto na resolução
        LIKWID_MARKER_STOP("GaussSeidel"); // Marca o fim da região medida

        // Imprime solução
        for (int i = 0; i < n; ++i)
            printf(" %23.15e", Y[i]);
        printf("\n");
        printf("%d\n", numIt);
        printf(" %23.15e\n", normaL2);
        printf(" %16.8e\n", t1);

        // Libera memória alocada para o sistema e solução
        free(Y);
        free(sl->D); 
        free(sl->Di); 
        free(sl->Ds); 
        free(sl->B); 
        free(sl);
    }

    //Finaliza os marcadores do LIKWID
    LIKWID_MARKER_CLOSE;
    return 0;
}