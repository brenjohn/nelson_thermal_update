#include "microphysics_dabundances.h"

/* The derivatives for n cells. All arrays are row-major: X is n x N_VARS (u, T), params n x N_PARAMS and dF
   n x N_VARS x N_DERIVATIVES. */
void batch_dabundances(int n, const double *X, const double *params, double *dF) {
    for (int i = 0; i < n; i++)
        microphysics_dabundances((const SolveVars *)(X + i * N_VARS), (const Params *)(params + i * N_PARAMS),
                                 (double (*)[N_DERIVATIVES])(dF + i * N_VARS * N_DERIVATIVES));
}
