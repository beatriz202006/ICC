#include <stdio.h>
#include <stdlib.h>
#include <fenv.h>
#include "edo.h"
#include "gauss_seidel.h"
#include "utils.h"

int main() {
    // Força arredondamento para baixo
    fesetround(FE_DOWNWARD);

    int n;
    real_t a, b, ya, yb, p, q;
    EDo edo;
    real_t r1, r2, r3, r4;

    // Lê dados iniciais
    if (scanf("%d", &n) != 1) return 1;
    if (scanf("%lf %lf", &a, &b) != 2) return 1;
    if (scanf("%lf %lf", &ya, &yb) != 2) return 1;
    if (scanf("%lf %lf", &p, &q) != 2) return 1;

    edo.n = n;
    edo.a = a;
    edo.b = b;
    edo.ya = ya;
    edo.yb = yb;
    edo.p = p;
    edo.q = q;

    // Para cada linha de coeficientes r1 r2 r3 r4
    while (scanf("%lf %lf %lf %lf", &r1, &r2, &r3, &r4) == 4) {
        edo.r1 = r1;
        edo.r2 = r2;
        edo.r3 = r3;
        edo.r4 = r4;

        Tridiag *sl = genTridiag(&edo);

        prnEDOsl(&edo);

        real_t *Y = (real_t *)calloc(n, sizeof(real_t));
        int maxIt = 100;
        real_t tol = 1e-5;
        real_t normaL2;
        int numIt;

        rtime_t t0 = timestamp();
        gaussSeidel(sl, Y, maxIt, tol, &normaL2, &numIt);
        rtime_t t1 = timestamp() - t0;

        // Imprime solução
        for (int i = 0; i < n; ++i)
            printf(" %23.15e", Y[i]);
        printf("\n");
        printf("%d\n", numIt);
        printf(" %23.15e\n", normaL2);
        printf(" %16.8e\n", t1);

        free(Y);
        free(sl->D); free(sl->Di); free(sl->Ds); free(sl->B); free(sl);
    }
    return 0;
}