#pragma once

#define JACO_HAS_VAR_u
#define JACO_HAS_VAR_T
#define JACO_HAS_PARAM_A_V
#define JACO_HAS_PARAM_Delta_t
#define JACO_HAS_PARAM_G_0
#define JACO_HAS_PARAM_dv_dr
#define JACO_HAS_PARAM_n_Htot
#define JACO_HAS_PARAM_pdv_work
#define JACO_HAS_PARAM_tau_uv
#define JACO_HAS_PARAM_u_initial
#define JACO_HAS_PARAM_x_CH
#define JACO_HAS_PARAM_x_CO
#define JACO_HAS_PARAM_x_C_tot
#define JACO_HAS_PARAM_x_Cplus
#define JACO_HAS_PARAM_x_H3plus
#define JACO_HAS_PARAM_x_HCOplus
#define JACO_HAS_PARAM_x_H_2
#define JACO_HAS_PARAM_x_Heplus
#define JACO_HAS_PARAM_x_Mg_tot
#define JACO_HAS_PARAM_x_Mgplus
#define JACO_HAS_PARAM_x_OH
#define JACO_HAS_PARAM_x_O_tot
#define JACO_HAS_PARAM_y
#define JACO_PARAM_NAMES {"A_V", "Delta_t", "G_0", "dv_dr", "n_Htot", "pdv_work", "tau_uv", "u_initial", "x_CH", "x_CO", "x_C_tot", "x_Cplus", "x_H3plus", "x_HCOplus", "x_H_2", "x_Heplus", "x_Mg_tot", "x_Mgplus", "x_OH", "x_O_tot", "y"}

enum SolveVarIndex { IDX_u = 0, IDX_T = 1, N_VARS = 2 };
enum ParamIndex { PARAM_A_V = 0, PARAM_Delta_t = 1, PARAM_G_0 = 2, PARAM_dv_dr = 3, PARAM_n_Htot = 4, PARAM_pdv_work = 5, PARAM_tau_uv = 6, PARAM_u_initial = 7, PARAM_x_CH = 8, PARAM_x_CO = 9, PARAM_x_C_tot = 10, PARAM_x_Cplus = 11, PARAM_x_H3plus = 12, PARAM_x_HCOplus = 13, PARAM_x_H_2 = 14, PARAM_x_Heplus = 15, PARAM_x_Mg_tot = 16, PARAM_x_Mgplus = 17, PARAM_x_OH = 18, PARAM_x_O_tot = 19, PARAM_y = 20, N_PARAMS = 21 };

typedef union {
    struct {
        double u;
        double T;
    };
    double data[2];
} SolveVars;

typedef union {
    struct {
        double A_V;  /* visual extinction of the photoreactions and the dust temperature [mag] */
        double Delta_t;  /* time step [s] */
        double G_0;  /* FUV field relative to the Habing (1968) field [Habing] */
        double dv_dr;  /* line-width (velocity) gradient of the C I escape probabilities, NL99: 1 km/s/pc [s^-1] */
        double n_Htot;  /* number density of H nuclei [cm^-3] */
        double pdv_work;  /* compressional (PdV) heating rate [erg cm^-3 s^-1] */
        double tau_uv;  /* UV optical depth of the photoelectric and surface terms */
        double u_initial;  /* specific internal energy at the start of the step [erg g^-1] */
        double x_CH;  /* CH per H nucleus */
        double x_CO;  /* CO per H nucleus */
        double x_C_tot;  /* C nuclei per H nucleus, all species */
        double x_Cplus;  /* C+ per H nucleus */
        double x_H3plus;  /* H3+ per H nucleus */
        double x_HCOplus;  /* HCO+ per H nucleus */
        double x_H_2;  /* H_2 per H nucleus */
        double x_Heplus;  /* He+ per H nucleus */
        double x_Mg_tot;  /* Mg nuclei per H nucleus, all species */
        double x_Mgplus;  /* Mg+ per H nucleus */
        double x_OH;  /* OH per H nucleus */
        double x_O_tot;  /* O nuclei per H nucleus, all species */
        double y;  /* He nuclei per H nucleus */
    };
    double data[21];
} Params;

void microphysics_func_jac(const SolveVars *vars, const Params *params, SolveVars *rhs, double jac[N_VARS][N_VARS]);

/* Outputs: quantities microphysics_outputs() evaluates at a state, without derivatives. At the converged
   state Delta_t times each is its integral over the backward-Euler step, so the host applies it directly. */
enum OutputIndex { N_OUTPUTS = 0 };
typedef union { double data[1]; } Outputs; /* the model has no outputs */
void microphysics_outputs(const SolveVars *vars, const Params *params, Outputs *out);

/* Solver metadata. Per solve variable, in IDX_ order: backward-Euler term (1) or steady state (0), the
   floor, ceiling and scale of a species' abundance per H nucleus (u and T take their limits from the
   solver settings), charge. */
#define JACO_HAS_SOLVER_METADATA
#define JACO_VAR_TIME_DEPENDENT_T
#define JACO_VAR_TIME_DEPENDENT_INIT {0, 1}
#define JACO_VAR_FLOOR_INIT {0.0, 0.0}
#define JACO_VAR_CEILING_INIT {1e+300, 1e+300}
#define JACO_VAR_SCALE_INIT {1.0, 1.0}
#define JACO_VAR_CHARGE_INIT {0, 0}
/* kind of each variable: JACO_KIND_ENERGY (u), _GAS_TEMPERATURE (T), _SPECIES (a material or trace abundance), _RADIATION (a band, per H nucleus), _TEMPERATURE (determined by the steady state of its row, e.g. the dust's) */
#define JACO_KIND_ENERGY 0
#define JACO_KIND_GAS_TEMPERATURE 1
#define JACO_KIND_SPECIES 2
#define JACO_KIND_RADIATION 3
#define JACO_KIND_TEMPERATURE 4
#define JACO_VAR_KIND_INIT {JACO_KIND_ENERGY, JACO_KIND_GAS_TEMPERATURE}
/* PARAM_ index of each variable's start-of-step value (-1: none; T's is u_initial) */
#define JACO_VAR_INITIAL_PARAM_INIT {-1, PARAM_u_initial}
/* Time-dependent species: {IDX_ of the species, PARAM_ of its start-of-step value} */
#define JACO_N_TD_SPECIES 0
#define JACO_TD_SPECIES_INIT { \
    {-1, -1} /* none */, \
}
/* Budgets of the eliminated abundances that are affine in the solve variables:
   total - sum_t w[t] x[k[t]] >= 0, as {total, total PARAM_ index or -1 for the constant, nterm, {IDX_...},
   {w...}} */
#define JACO_N_BUDGETS 0
#define JACO_BUDGET_MAX_TERMS 1
#define JACO_BUDGETS_INIT { \
    {0.0, -1, 0, {0}, {0.0}} /* none */, \
}
