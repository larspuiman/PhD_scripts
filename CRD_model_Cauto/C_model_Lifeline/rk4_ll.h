
double *rk4vec(double t0, double u0[12], double dt, double conc_ref[12], double rates_IC_ref[11], double elasticity_matrix[11][12], double Np);
double *conc_derivatives(double conc_act[12], double conc_ref[12], double rates_IC_ref[11], double elasticity_matrix[11][12], double rate_vec[], double Np);