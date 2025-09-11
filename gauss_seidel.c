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
    for (it = 0; it <= maxIt; it++) {

        //Copia X para X_old:
        for (int i = 0; i < n; i++)
            X_old[i] = X[i];
        
        //Iteração de Gauss-Seidel:
        for (int i = 0; i < n; i++) {
            real_t soma = sl->B[i];

            //Primeira equacao nao tem diagonal inferior
            if (i > 0)
                soma -= sl->Di[i] * X[i-1];
            
            //Ultima equacao nao tem diagonal superior
            if (i < n-1)
                soma -= sl->Ds[i] * X[i+1];
            
            X[i] = soma / sl->D[i];
        }

        //Calcula norma L2 do erro:
        *norma2 = 0.0;
        for (int i = 0; i < n; i++)
            *norma2 += (X[i] - X_old[i]) * (X[i] - X_old[i]);
        
        *norma2 = sqrt(*norma2);

        //Se convergiu: 
        if (*norma2 <= tol) 
            break;
        
    }

    *numIt = it;
    free(X_old);
    return 0; //Convergiu com sucesso
}