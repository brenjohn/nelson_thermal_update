#include "microphysics_func_jac.h"

/* The generated function for n cells. All arrays are row-major: X is n x N_VARS (u, T), params n x N_PARAMS,
   rhs n x N_VARS and jac n x N_VARS x N_VARS. */
void batch_func_jac(int n, const double *X, const double *params, double *rhs, double *jac) {
    for (int i = 0; i < n; i++)
        microphysics_func_jac((const SolveVars *)(X + i * N_VARS), (const Params *)(params + i * N_PARAMS),
                              (SolveVars *)(rhs + i * N_VARS), (double (*)[N_VARS])(jac + i * N_VARS * N_VARS));
}
