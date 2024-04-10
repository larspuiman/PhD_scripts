#include "udf.h"
#include "dpm.h"
#include "mem.h"
#include "math.h"
#include "stdio.h"

#if RP_DOUBLE == 0
#error "UDF requires double precision mode"
#endif

/* define the number of compounds and reactions (equals the size of eps, rates) */
int num_compounds = 12;
int num_rxns = 11;
float eps = 1e-15;

/* Define global constants */
real MW_CO = 28.01;		/* Molecular mass: g/mol */
real MW_H2 = 2.016;		/* Molecular mass: g/mol */
real MW_CO2 = 44.01;	/* Molecular mass: g/mol */

/* Functions */
int thread_counter();
void get_rates_matrix_fun();
get_minimum_conc();
double * read_matrix(int rows, int cols, const char* filename);
double * conc_derivatives(double conc_act[num_compounds],  double dt, double NP_local, double* dc_dt, double rates_IC_ref[num_rxns], double conc_ref[num_compounds], double conc_min[num_compounds], double elasticity_matrix[num_rxns][num_compounds]);
double * read_rates(int n_vals, const char* filename);
double * read_conc(int n_vals, const char* filename);
double * rk4vec ( double t0, double* u0, double dt, double NP_local, double rates_IC_ref[num_rxns], double conc_ref[num_compounds], double conc_min[num_compounds], double elasticity_matrix[num_rxns][num_compounds]);
double * rates_fun(double conc_act[num_compounds], double* rates, double rates_IC_ref[num_rxns], double conc_ref[num_compounds], double conc_min[num_compounds], double elasticity_matrix[num_rxns][num_compounds] );
