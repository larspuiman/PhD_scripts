
void rk4vec(double t0, double u0[12], double dt, double conc_ref[12], double rates_IC_ref[11], double elasticity_matrix[11][12], double u[12], double CX);
void conc_derivatives(double conc_act[12], double conc_ref[12], double rates_IC_ref[11], double elasticity_matrix[11][12], double dc_dt_new[12], double CX);