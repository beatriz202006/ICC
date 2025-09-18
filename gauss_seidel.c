/*Beatriz Pontes Camargo 
GRR 20242966
*/
#include <math.h>
#include <stdlib.h>
#include "gauss_seidel.h"
#include "utils.h"

//Resolve SL tridiagonal pelo método de Gauss-Seidel
int gaussSeidel(Tridiag *sl, real_t *X, int maxIt, real_t tol, real_t *norma2, int *numIt)
{
    int n = sl->n;
    real_t *X_old = (real_t *)malloc(n * sizeof(real_t));

    //Verifica alocação
    if (!X_old) 
        return 1;
    
    int it;
    real_t norma_res = 0.0;

    for (it = 0; it <= maxIt; it++) {

        //Copia X para X_old:
        for (int i = 0; i < n; i++)
            X_old[i] = X[i];
        
        //Iteração de Gauss-Seidel:
        for (int i = 0; i < n; i++) {
            real_t soma = sl->B[i];

            //Primeira equação não tem diagonal inferior
            if (i > 0)
                soma -= sl->Di[i] * X[i-1];
            
            //Última equação não tem diagonal superior
            if (i < n-1)
                soma -= sl->Ds[i] * X[i+1];
            
            X[i] = soma / sl->D[i];
        }

        norma_res = 0.0;
        //Critério de parada: normaL2 do resíduo => calcula o resíduo a cada iteração
        for (int i = 0; i < n; i++) {
            real_t Ax = sl->D[i] * X[i];
            if (i > 0)
                Ax += sl->Di[i] * X[i-1];
            if (i < n-1)
                Ax += sl->Ds[i] * X[i+1];
            real_t Ri = sl->B[i] - Ax;
            norma_res += Ri * Ri;
        }
        norma_res = sqrt(norma_res);

        // Critério de parada: norma L2 do resíduo
        if (norma_res <= tol)
            break;
            
        }

    *numIt = it+1;
    *norma2 = norma_res;

    free(X_old);
    return 0; //Convergiu com sucesso
}