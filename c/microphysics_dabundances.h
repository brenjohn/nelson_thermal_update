#pragma once
#include "microphysics_func_jac.h"

/* Derivatives of the residuals of microphysics_func_jac with respect to the abundances (per H nucleus) of the species
   the model takes as inputs, and to u_initial, at the solve variables in vars. The columns are in the order of
   DerivativeIndex. */
enum DerivativeIndex { D_x_H_2 = 0, D_x_H3plus = 1, D_x_Heplus = 2, D_x_CH = 3, D_x_OH = 4, D_x_CO = 5, D_x_Cplus = 6, D_x_HCOplus = 7, D_x_Mgplus = 8, D_u_initial = 9, N_DERIVATIVES = 10 };

void microphysics_dabundances(const SolveVars *vars, const Params *params, double dF[N_VARS][N_DERIVATIVES]);
