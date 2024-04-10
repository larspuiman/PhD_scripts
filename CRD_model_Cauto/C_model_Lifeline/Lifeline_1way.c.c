#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "rk4_ll.h"


/* FUNCTION to read a certain number of rows and columns from a matrix in a .dat file, specified by filename*/
double *read_matrix(int rows, int cols, const char* filename)
{
    double matrix[rows][cols];
    double *a;

    a = ( double * ) malloc ( rows * cols * sizeof ( double ) ); 
    FILE *file_pointer;
    file_pointer = fopen (filename, "r");
    if (file_pointer == NULL)
        return 0;

    for(int i = 0; i < rows; ++i)
    {
        for(int j = 0; j < cols; ++j)
            fscanf(file_pointer, "%lf", &matrix[i][j]);
    }

    for(int i = 0; i < rows; ++i)
    {
        for(int j = 0; j < cols; ++j)
            {
                int offset = i *cols + j;
                a[offset] = matrix[i][j];
            }
    }

    fclose (file_pointer); 
    return a; 
}

/* FUNCTION to write a matrix */
double writematrix(int rows, int cols, double (*matrix)[cols], const char* filename)
{
    FILE *file_pointer;
    file_pointer = fopen (filename, "w+");
    if (file_pointer == NULL)
        return 0;

    for(int i = 0; i < rows; ++i)
    {
        for(int j = 0; j < cols; ++j)
            fprintf(file_pointer, "%e ", matrix[i][j]);
        fprintf(file_pointer, "\n");
    }

    fclose (file_pointer); 
    return 1; 
}

/* FUNCTION to read a vector of length n_vals from a .dat file, specified by filename*/
double *read_conc(int num_compounds, const char* filename)
{
    double conc[num_compounds];
    double *a;
    a = ( double * ) malloc ( num_compounds * sizeof ( double ) );
    FILE *file_pointer;
    file_pointer = fopen (filename, "r");
    if (file_pointer == NULL)
        return 0;

    for(int i = 0; i < num_compounds; ++i)
    {
        fscanf(file_pointer, "%lf", &conc[i]);
        a[i] = conc[i];
    }

    fclose (file_pointer); 
    return a; 
}

/* FUNCTION to read a vector of length n_vals from a .dat file, specified by filename*/
double *read_rates(int num_rxns, const char* filename)
{
    double rates[num_rxns];
    double *a;
    a = ( double * ) malloc ( num_rxns * sizeof ( double ) );
    FILE *file_pointer;
    file_pointer = fopen (filename, "r");
    if (file_pointer == NULL)
        return 0;

    for(int i = 0; i < num_rxns; ++i)
    {
        fscanf(file_pointer, "%lf", &rates[i]);
        a[i] = rates[i];
    }

    fclose (file_pointer); 
    return a; 
}

/* FUNCTION to read a vector of length n_vals from a .dat file, specified by filename*/
double *read_lifeline(int n_vals, const char* filename)
{
    double vector[n_vals];
    double *a;
    a = ( double * ) malloc ( n_vals * sizeof ( double ) );
    FILE *file_pointer;
    file_pointer = fopen (filename, "r");
    if (file_pointer == NULL)
        return 0;

    for(int i = 0; i < n_vals; ++i)
    {
        fscanf(file_pointer, "%lf", &vector[i]);
        a[i] = vector[i];
    }

    fclose (file_pointer); 
    return a; 
}


/* FUNCTION to write a vector */
double writevector(int n_vals, double (*vector), const char* filename)
{
    FILE *file_pointer;
    file_pointer = fopen (filename, "w+");
    if (file_pointer == NULL)
        return 0;

    for(int i = 0; i < n_vals; ++i)
        fprintf(file_pointer, "%e \n ", vector[i]);
    
    fclose (file_pointer); 
    return 1; 
}

double *rk4vec ( double t0, double *u0, double dt, double conc_ref[12], double rates_IC_ref[11], double elasticity_matrix[11][12], double Np )

/******************************************************************************/
/*
  Purpose:
 
    RK4VEC takes one Runge-Kutta step for a vector ODE.

  Discussion:

    It is assumed that an initial value problem, of the form

      du/dt = f ( t, u )
      u(t0) = u0

    is being solved.

    If the user can supply current values of t, u, a stepsize dt, and a
    function to evaluate the derivative, this function can compute the
    fourth-order Runge Kutta estimate to the solution at time t+dt.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, double T0, the current time.

    Input, int M, the spatial dimension.

    Input, double U0[M], the solution estimate at the current time.

    Input, double DT, the time step.

    Input, double *F ( double T, int M, double U[] ), a function which evaluates
    the derivative, or right hand side of the problem.

    Output, double RK4VEC[M], the fourth-order Runge-Kutta solution estimate
    at time T0+DT.
*/
{
  double *f0;
  double *f1;
  double *f2;
  double *f3;
  int i;
  double t1;
  double t2;
  double t3;
  double *u;
  double *u1;
  double *u2;
  double *u3;

  int num_compounds = 12;
  double rate_vec[num_compounds];
/*
  Get four sample values of the derivative.
*/
  f0 = conc_derivatives ( u0 , conc_ref, rates_IC_ref, elasticity_matrix, rate_vec, Np ) ;

  t1 = t0 + dt / 2.0;
  u1 = ( double * ) malloc ( num_compounds * sizeof ( double ) );
  for ( i = 0; i < num_compounds ; ++i )
  {
    u1[i] = u0[i] + dt * f0[i] / 2.0;
  }
  f1 = conc_derivatives (  u1,  conc_ref, rates_IC_ref, elasticity_matrix, rate_vec, Np );

  t2 = t0 + dt / 2.0;
  u2 = ( double * ) malloc ( num_compounds * sizeof ( double ) );
  for ( i = 0; i < num_compounds; ++i )
  {
    u2[i] = u0[i] + dt * f1[i] / 2.0;
  }
  f2 = conc_derivatives ( u2,  conc_ref, rates_IC_ref, elasticity_matrix, rate_vec, Np );

  t3 = t0 + dt;
  u3 = ( double * ) malloc ( num_compounds * sizeof ( double ) );
  for ( i = 0; i < num_compounds; ++i )
  {
     u3[i] = u0[i] + dt * f2[i];
  }
  f3 = conc_derivatives ( u3, conc_ref, rates_IC_ref, elasticity_matrix, rate_vec, Np );
/*
  Combine them to estimate the solution.
*/
  u = ( double * ) malloc ( num_compounds * sizeof ( double ) );
  double eps = 1e-15;
//   double min_concv[] =  { [0] = 1e-10, [1] = 1e-10, [2 ... 5] = 0.1, [6] = 1e-9, [7] = 1e-12, [8] = 1e-12, [9] = 1e-12, [10] = 1e-12, [11] = 1e-15};
  double min_concv[] =  { [0] = 1e-10, [1] = 1e-10, [2 ... 5] = 0.1, [6] = 1e-8, [7] = 1e-12, [8] = 1e-12, [9] = 1e-12, [10] = 1e-4, [11] = 1e-12};
//   double min_concv[] =  { [0] = 1e-5, [1] = 1e-5, [2 ... 5] = 0.1, [6] = 1e-3, [7] = 1e-3, [8] = 1e-3, [9] = 1e-9, [10] = 1e-8, [11] = 1e-9};

  /* calculate new conccentrations */
  for ( i = 0; i < num_compounds; ++i )
  {
    u[i] = u0[i] + dt * ( f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i] ) / 6.0;

    /* in case the concentration is lower than the minimum allowable one, give the conc that value */
    if ( (u[i] < min_concv[i]) )
        {u[i] = min_concv[i]*0.9;}
  }

/*
  Free memory.
*/
  free ( f0 );
  free ( f1 );
  free ( f2 );
  free ( f3 );
  free ( u1 );
  free ( u2 );
  free ( u3 );

  return u;
}

/* FUNCTION to calculate derivatives */
double *conc_derivatives(double conc_act[12], double conc_ref[12], double rates_IC_ref[11], double elasticity_matrix[11][12], double rate_vec[], double Np )
{    
    /* provide  number of compounds and reactions */
    int num_compounds = 12;                    
    int num_rxns =  11;       

    double *dc_dt = calloc(num_compounds, sizeof(double));
    double eps = 1e-8;
    // double min_concv[] = { [0] = 1e-10, [1] = 1e-10, [2 ... 5] = 0.1, [6] = 1e-9, [7] = 1e-12, [8] = 1e-12, [9] = 1e-12, [10] = 1e-12, [11] = 1e-15};
    double min_concv[] =  { [0] = 1e-10, [1] = 1e-10, [2 ... 5] = 0.1, [6] = 1e-8, [7] = 1e-12, [8] = 1e-12, [9] = 1e-12, [10] = 1e-4, [11] = 1e-12};
    /* definitions */
    double conc_log[num_compounds];
    double rates_ll_cor[num_rxns];
    double rates_IC[num_rxns];
    double r_HAcD, dcHdt;

    /* constants */
    double conc_eff[num_compounds];
    int neg_conc[12];
    double cell_surf_area = 321;                     /* surface area of cells [m2 / molx] */
    double k_HAcD = 3.85e-5;                         /* diffusivity of acetic acid [m / h] */
    double pH_EC = 5;                                /* extracelluar pH ~ variable in the CFD model */
    double pKa_HAc = 4.756;                          /* Acetic acid pKa (source: wikipedia) */
    double CX = 100;                                 /* biomass concentration LATER [molx/m3] */
    double Vx = 5.89e-5;                             /* Bacteria molar volume [m3/molx] */

    /* check whether the concentrations are too low (c < eps and store then) */
    for(int i = 0; i < num_compounds; ++i){
        if (conc_act[i] < min_concv[i])
        {conc_eff[i] = min_concv[i];}
        else
        {conc_eff[i] = conc_act[i];}
    }

    /* get the natural logarithm of the IC conc vs its ref [ - ] */
    for(int i = 0; i < num_compounds; ++i)
        conc_log[i] = log(conc_eff[i]/conc_ref[i]);

    /* calculate the linlog corrected rates [ - ] */
    for( int i = 0; i < num_rxns; ++i)
    { 
        double sum_elast = 0;
        for(int j = 0; j < num_compounds; ++j)
        {
            int sign = 1; /* multiply by -1 when the conc are too low so that the rates get zero in the end */
            // if ((elasticity_matrix[i][j] < 0) && (conc_log[j] < -10) && (i != 1) )
            //     {sign = -1;}
            sum_elast += elasticity_matrix[i][j] * sign * conc_log[j];
        }
        
        rates_ll_cor[i]= 1 + sum_elast;

        // /* the irreversible reactions (ACAS, AcS, Rnf, Ana, AcX cannot become negative)*/
        if (rates_ll_cor[i] < 0  && (i == 1 || i == 4 || i == 8 || i == 9 || i == 10 ) ) 
            {rates_ll_cor[i] = 0;}
    }

    /* calculate the total rates of intracellular reactions [mol/molx/h ]*/
    for(int i = 0; i < num_rxns; ++i)
    { 
        rates_IC[i]= rates_IC_ref[i]*rates_ll_cor[i];
    }

    /* calculate the rate of acetate back diffusion [mol/molx/h] */
    double c_HAc_EC = conc_act[3]*pow(10,pKa_HAc)/(pow(10,pH_EC)+pow(10,pKa_HAc)); /* conc of HAc (extracellular) [mol/m3l]*/
    r_HAcD = cell_surf_area*k_HAcD*c_HAc_EC;                                      /* back diffusion rate [mol/molx/h] */

    /* calculate the differentials */
    double mu = rates_IC[9];                                         /* biomass growth rate */
        /* for cIC = cEC [mol/m3_L/h] */
    dc_dt[0] = -1.*(rates_IC[1] + rates_IC[2])   * CX ;              /* CO uptake rate [mol/m3/h] (negative for consumption) */
    dc_dt[1] = -2.*(rates_IC[3])                 * CX ;              /* H2 uptake rate [mol/m3/h] (negative for consumption) */
    dc_dt[2] =    (rates_IC[2] - rates_IC[0])   * CX ;              /* CO2 uptake rate [mol/m3/h] (negative for consumption)*/
    dc_dt[3] =    (rates_IC[10] - r_HAcD)       * CX ;              /* Net acetate excretion rate [mol/m3/h] */
    dc_dt[4] =    (rates_IC[5])                 * CX ;              /* EtOH prod rate [mol/m3/h] */
    dc_dt[5] =    (rates_IC[6])                 * CX ;              /* BDO prod rate [mol/m3/h] */
        /* only compounds that remain in cell [mol/m3_cell/h] */
    dc_dt[6] =  (rates_IC[1] -   rates_IC[4] -    2.*rates_IC[6] - 1./2.*rates_IC[9])               * (1.0/Vx) - mu*conc_act[6];   /* AcCoA balance [mol/m3c/h] */
    dc_dt[7] =  (rates_IC[8] - 2.*rates_IC[1] -     rates_IC[5]  - 1./2.*rates_IC[6] - rates_IC[7]) * (1.0/Vx) - mu*conc_act[7];   /* NADH balance [mol/m3c/h] */
    dc_dt[8] =  (rates_IC[3] + 2.*rates_IC[7] - 1./2.*rates_IC[0]  -     rates_IC[1])               * (1.0/Vx) - mu*conc_act[8];   /* NADPH balance [mol/m3c/h] */
    dc_dt[9] =  (rates_IC[4] +   r_HAcD   -     rates_IC[5]  -     rates_IC[10])                 * (1.0/Vx) - mu*conc_act[9];   /* Ac_i_in balance [mol/m3c/h] */
    dc_dt[10] = (rates_IC[0] -   rates_IC[1])                                                    * (1.0/Vx) - mu*conc_act[10];  /* Formate balance [mol/m3c/h] */
    dc_dt[11] = (-0.5*rates_IC[0] + rates_IC[1] + rates_IC[2] + rates_IC[3] - 
                rates_IC[5] - rates_IC[6] - rates_IC[7] - rates_IC[8])                           * (1.0/Vx) - mu*conc_act[11];  /* Ferredoxin balance [mol/m3c/h] */

    dcHdt = -1*(-0.5*rates_IC[0] - 2*rates_IC[1] + 2*rates_IC[2] + 3*rates_IC[3] + 1*rates_IC[4] - 4*rates_IC[5] 
            - 5./2.*rates_IC[6] - 1*rates_IC[7 ] - 3*rates_IC[8] + 1./4.*rates_IC[9] + r_HAcD);                             /* H+ rate [mol/m3/h] */

    /* store the rates via the rate vector */
    for(int i = 0; i < num_compounds; ++i)
    { 
        if(i == num_rxns)
            {rate_vec[i] = r_HAcD;}
        else
            {rate_vec[i] = rates_IC[i];}
    }
    
    /* implement the constraints that the max dc/dt difference with the maxrate is governing the stoichiomtery*/
    double max_dcdt[num_compounds];                         /* the compound specific max dc/dt (mol/m3/h) */
    double dcdt_ratio[num_compounds];                       /* stoichiometries */
    double dcdt_diff[num_compounds];                        /* the difference between the calculated dcdt and its max (should not be negative!) */
    double dt = 1e-4/3600;                                  /* timestep, should be implemented via a function argumetn */
    int max_index = 0; double max_rate = 0;                 /* defs for maximum */
    int zero_rate = 0;                                      /* flag to store if some rates are negative and should be corrected therefore */
    for(int i = 0; i < num_compounds; ++i)
    {
        if (i < 6)
            {max_dcdt[i] = (conc_eff[i]-min_concv[i]*0.9)/(Np*dt);}        /* maximum rate EC compounds [mol/m3L/part/h] */
        else
            {max_dcdt[i] = (conc_eff[i]-min_concv[i]*0.9)/dt;}                   /* maximum rate IC compounds [mol/m3c/h] */
        dcdt_diff[i] = dc_dt[i] - max_dcdt[i];              /* diff between current rate and its max */

        /* determine the index for the compound with the largest rate difference */
        if (dcdt_diff[i] > max_rate)
            {
                max_rate = dcdt_diff[i];
                max_index = i;
            }
        if (dc_dt[i] > max_dcdt[i])                             /* if some rates are zero */
            {zero_rate = 1;}
    }

    if (zero_rate == 1)
        {
            /* make the dcdt with the max_index the new reference rate  */
            double dc_dt_ref = max_dcdt[max_index];

            /* correct the other dcdt with that factor, considering the current stoichiometry (dcdt-ratio's)  */
            for(int i = 0; i < num_compounds; ++i)
            {
                dc_dt[i] = dc_dt_ref/(dc_dt[max_index]) * dc_dt[i];               /* new dcdt */
            }
        }

    return dc_dt;
}

/* FUCTION to calculate averages over integration (Eulerian source term) */
double int_average(double eulerian_timestep, double storage_freq, double conc[])
{
    double time_per_storage = eulerian_timestep/storage_freq;
    double sum_conc_dt = 0;
    for(int i = 0; i < (storage_freq+1); i++)
    {
        sum_conc_dt += conc[i]*time_per_storage;
    }
    return sum_conc_dt/eulerian_timestep;
}

int main(void)
{
    /* define the number of compounds and reactions (equals the size of eps, rates) */
    int num_compounds = 12;
    int num_rxns = 11;
    int num_time_data = 13339;

    /* define the vectors for concentration, rates and elasticity matrix */
//    double conc_act[num_compounds];               /* Actual intracellular conc */
    double conc_av[num_compounds];                /* Average concentration during Lagrangian timestep */
    double elasticity_matrix[num_rxns][num_compounds];
    double conc_start[num_compounds];

    double *conc_ref = calloc(num_compounds, sizeof(double));
    double *conc_ini = calloc(num_compounds, sizeof(double));
    double *rates_IC_ref = calloc(num_compounds, sizeof(double));
    double *matrix_temp = calloc(num_compounds*num_rxns, sizeof(double));
    
    double *CO_lifeline = calloc(num_time_data, sizeof(double));
    double *H2_lifeline = calloc(num_time_data, sizeof(double));
    double *CO2_lifeline = calloc(num_time_data, sizeof(double));
    double *t_lifeline = calloc(num_time_data, sizeof(double));
    double *Np_lifeline = calloc(num_time_data, sizeof(double));


    /* stuff that gets updated all the time */
    double *dc_dt_ini = calloc(num_compounds, sizeof(double));
    double *conc_act = calloc(num_compounds, sizeof(double)); 

    double dcdt_av[num_compounds];                /* Average dc/dt during Lagrangian timestep */
    double rate_vec[num_compounds];               /* Rates during Lagrangian timestep */
    
    /* define clock variables */
    clock_t start_t, end_t; 
    double total_t;

    /* read the concentrations of the reference state [ mol/m3 ]  */
    conc_ref = read_conc(num_compounds, "Ref_conc.dat");

    /* Get the initial concentrations  [ mol/m3 ]  */
    conc_ini = read_conc(num_compounds, "Ini_conc.dat");

    /* read the rates of the reference state [ mol/molx/h ] */
    rates_IC_ref = read_rates(num_rxns, "Ref_rates.dat");

    /* Read the lifeline  [ mol/m3 ]  */
    CO_lifeline =  read_lifeline(num_time_data, "CO_lifeline.dat");
    H2_lifeline = read_lifeline(num_time_data, "H2_lifeline.dat");
    CO2_lifeline = read_lifeline(num_time_data, "CO2_lifeline.dat");
    t_lifeline = read_lifeline(num_time_data, "t_lifeline.dat");
    Np_lifeline = read_lifeline(num_time_data, "Np_lifeline.dat");

    /* read the matrix */
    matrix_temp = read_matrix(num_rxns, num_compounds, "eps_matrix.dat");
    for(int i = 0; i < num_rxns; ++i)
    {
        for(int j = 0; j < num_compounds; ++j)
            {
                int offset = i * num_compounds + j;
                elasticity_matrix[i][j] = matrix_temp[offset];
            }
    }    


    /* Variables for the integration */
    double lagrangian_timestep = 1e-4/3600;                                            /* timestep in [h] */
    double int_time = 1000.0/3600;                                              /* timestep in [h] */
    int num_timesteps = (int) (int_time/lagrangian_timestep);                   /* timesteps in the Lagrangian integration */
    double eulerian_timestep = 0.075/3600;                                                  /* eulerian timestep */
    int timestep_diff = (int) (eulerian_timestep/lagrangian_timestep)+1; 
    int storage_freq = 1000;
    int storage_points = num_timesteps/storage_freq;
  
    /* set CO, H2, CO2 as initial conc from lifeline */
    conc_ini[0] = CO_lifeline[0];
    conc_ini[1] = H2_lifeline[0];
    conc_ini[2] = CO2_lifeline[0];

    double *dc_dt_new = calloc(num_compounds, sizeof(double)); 
    double *conc_new = calloc(num_compounds, sizeof(double)); 
    double *dc_dt_old = calloc(num_compounds, sizeof(double));

    double conc_ini_etoh_1 = conc_ini[4];
    double conc_ref_etoh_1 = conc_ref[4];
    double rate_acx_ini = rates_IC_ref[10];
    double Np = Np_lifeline[0];

    /* prepare for running the lifeline */
    for (double c_Ac = 90; c_Ac <= 91; c_Ac += 5) 
    {
        // double c_Ac_switch = 39.0;
        // if (c_Ac < c_Ac_switch)
        //     {
        //     conc_ref[3] = 0.5;
        //     conc_ini[4] = 0.5;
        //     conc_ref[4] = 0.5;
        //     conc_ref[9] = 0.5;
        //     }
        // else
        //     {
            // conc_ref[3] = 60;
            // conc_ini[4] = 434.12;                  /* 20 g/L */
            // conc_ref[4] = conc_ref_etoh_1*434.12/119.7;
            
            // }
       // rates_IC_ref[10] = rate_acx_ini*1.05;    

        /* Allocate memory for storage */
        int store_index = 1;
        double current_time = 0;                                             /* maybe find something to increase resolution of float */
        double conc_store[num_compounds][storage_points+1];                     /* storage for the concentrations during integration */
        double dc_dt_store[num_compounds][storage_points+1];                    /* storage for the dcdt during integration */
        double time_data[storage_points+1];                                     /* vector of time data for export */ 
        start_t = clock();  
        
        

        double pH_EC = 5;                                /* extracelluar pH ~ variable in the CFD model */
        double pKa_HAc = 4.756;                          /* Acetic acid pKa (source: wikipedia) */

        double c_ac_d_out = c_Ac * pow(10, pH_EC - pKa_HAc) / (1 + pow(10, pH_EC - pKa_HAc) );
        double log_c_aco = log10(c_Ac - c_ac_d_out);
        double c_ac_d_in = pow(10, log_c_aco - (0. - 1.)*96.5*-0.0974/(8.314*310*log(10)) );
        // conc_ini[9] = c_ac_d_in;
        /* Calculate the derivatives using the derivate function */
        dc_dt_ini = conc_derivatives(conc_ini, conc_ref, rates_IC_ref, elasticity_matrix, rate_vec, Np );

        /* store t = 0 data */
        for(int i = 0; i < num_compounds; ++i)
        { 
            dc_dt_store[i][0] = dc_dt_ini[i];
            time_data[0] = current_time*3600;
            conc_act[i] = conc_ini[i];
            conc_start[i] = conc_act[i];
            if (i==3) {conc_act[i] = c_Ac;}
            conc_store[i][0] = conc_act[i];
        }
        /* set Ac concentration to be the iterator */
        conc_act[3] = c_Ac;

        printf("conc AcT = %f, conc Ac_IC = %f, rate AcX ref = %f\n", conc_act[3], conc_act[9], rates_IC_ref[10]);

        printf("Start loop - %d timesteps \n", num_timesteps);
        
        int lifeline_index = 0; int euler_counter = 0;
        for(int j = 0; j < num_timesteps; ++j) 
        {
                /* set CO, H2, CO2 as initial conc from lifeline, but take the next timestep into account */
                double t_proxy = (double) euler_counter/timestep_diff; 
                conc_act[0] = CO_lifeline[lifeline_index] * (1-t_proxy) +  CO_lifeline[lifeline_index+1] * t_proxy;
                conc_act[1] = H2_lifeline[lifeline_index] * (1-t_proxy) +  H2_lifeline[lifeline_index+1] * t_proxy;
                conc_act[2] = CO2_lifeline[lifeline_index] * (1-t_proxy) +  CO2_lifeline[lifeline_index+1] * t_proxy;
                Np = Np_lifeline[lifeline_index] * (1-t_proxy) +  Np_lifeline[lifeline_index+1] * t_proxy;

            /* time-index for lifeline */
            if (euler_counter == timestep_diff) {       /* we reached a new eulerian timestep - update concentrations  */
                lifeline_index += 1 ;
                euler_counter = 0; 
                // if (lifeline_index == 1334)
                //     {lifeline_index = 1;}

                for(int i = 0; i < num_compounds; ++i)
                {
                    conc_start[i] = conc_act[i];
                }

            }

            /* Calculate concentrations for current timestep */
            conc_new = rk4vec ( current_time, conc_act, lagrangian_timestep, conc_ref, rates_IC_ref, elasticity_matrix, Np);

            /* calculate the derivatives one more time for checking the results. The dcdt's are being used to check 
            whether the concentrations are in a good range - to prevent strange drops or so
                        the high and low concentration margin might be 100 times varying from the forward-euler solution */                            
            double* dc_dt_act =  (double*) calloc(num_compounds, sizeof(double)); 
            dc_dt_act = conc_derivatives(conc_act, conc_ref, rates_IC_ref, elasticity_matrix, rate_vec, Np);     /* [mol/m3/h] */
            double margin_factor = 100.0;                                       

            /* now lets do some temporary labelling based on the individual new concetration result */
            int nan_flag_count = 0; int good_result_flag_count = 0; int stiff_result_flag_count = 0; 
            for( int i = 0; i < num_compounds; ++i)
            {
                double high_margin = conc_act[i] + fabs(dc_dt_act[i]) * margin_factor * lagrangian_timestep/3600;
                double low_margin = conc_act[i] - fabs(dc_dt_act[i]) * margin_factor * lagrangian_timestep/3600;

                /* result is -nan or nan */
                if ( fabs(conc_new[i]) != fabs(conc_new[i]) )
                {
                    nan_flag_count += 1;
                }
                /* result is good - within the range predicted by forward Euler */
                else if (conc_new[i] < high_margin && conc_new[i] > low_margin )
                {
                    good_result_flag_count += 1;
                }   
                /* result is stiff or outside the range predicted by forward Euler */
                else if ( (conc_new[i] > high_margin || conc_new[i] < low_margin ) && conc_new[i] < 1e-5)
                {
                    stiff_result_flag_count += 1;
                }   

            }
               
            /* now the labelling of the general result. 
            Watch out, once flagged during the particle timestep, forever 
            flagged. So it might be a big overestimation */
            if (nan_flag_count >= 1)
            {
                for( int i = 0; i < num_compounds; ++i)
                {
                    conc_act[i] = conc_start[i];
                }
            }
            else if (good_result_flag_count > 10)
            {
                for( int i = 0; i < num_compounds; ++i)
                {
                    conc_act[i] = conc_new[i];
                }
            }

            
            free(dc_dt_act);




            for(int i = 0; i < num_compounds; ++i)
            { 
                /* data storage every X timesteps */
                if ((j+1)%(num_timesteps/storage_points) == 0) {
                    conc_store[i][store_index] = conc_new[i];
                    dc_dt_store[i][store_index] = dc_dt_old[i];
                }
                /* update the concentrations */
                conc_act[i] = conc_new[i];
            }
            /* update the current concentrations */

            /* calculate the new derivatives */
            dc_dt_new = conc_derivatives(conc_act, conc_ref, rates_IC_ref, elasticity_matrix, rate_vec, Np);
            for(int i = 0; i < num_compounds; ++i)
                dc_dt_old[i] = dc_dt_new[i];
            
            /* Update to next lagrangian timestep*/
            current_time += lagrangian_timestep;
            if ((j+1)%(num_timesteps/storage_points) == 0)
            {
                time_data[store_index] = current_time*3600;
                store_index += 1;
            }
            // free(conc_new); free(dc_dt_new); free(dc_dt_old);
            euler_counter += 1;
            // printf("end of iter, time = %f \n", time_data[j]);
        }    
        end_t = clock();    
        total_t =  (double)(end_t - start_t) / CLOCKS_PER_SEC;
        printf("Total time taken by CPU: %.12e\n", total_t  );

    
        /* write conc_data */
        char Filename_conc[150];
        char Filename_dcdt[150];
        char Filename_time[150];
        sprintf(Filename_conc, "Results/conc_data_1way_ll8.dat", (int) c_Ac);
        sprintf(Filename_dcdt, "Results/dcdt_data_1way_ll8.dat", (int) c_Ac);
        sprintf(Filename_time, "Results/t_data_1way_ll8.dat", (int) c_Ac);

        writematrix(num_compounds, storage_points+1, conc_store, Filename_conc);
        writematrix(num_compounds, storage_points+1, dc_dt_store, Filename_dcdt);
        writevector(storage_points+1, time_data, Filename_time);


    }




    



    return 0;
}

