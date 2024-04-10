
#include "udf.h"
#include "dpm.h"
#include "mem.h"
#include "math.h"
#include "stdio.h"
#include "unistd.h"
#include "dpm_types.h"
#include "prf.h" 

int i = 0;
int j = 0;

#define ZONE_ID 422

/* Define the number of compounds and reactions (equals the size of eps, rates) */
#define num_compounds 12
#define num_rxns 11                                 /* lazyness caused that sometimes rates + backdifussion has been defined as num_compounds... */                            

/* Do we want to write particle data every particle timestep? (useful for debugging but takes a lot of data) */
#define write_per_timestep 0    
#define write_per_while_timestep 0    
 
/* How many particles shall we track and write for? */
#define num_particles 80000
#define num_particles_write 20
#define num_particles_write_big 4000    
#define frequency_write_big 15                       /* write when   mod (timer/frequency) == 0   */
#define frequency_write_lifelines 15                 /* write when   mod (timer/frequency) == 0   */

/* Define global constants */
#define MW_CO 28.01   	                            /* Molecular mass [g/mol] */
#define MW_H2 2.016		                            /* Molecular mass [g/mol] */
#define MW_CO2 44.01	                            /* Molecular mass [g/mol] */
#define CX 150.0                                    /* biomass concentration  [molx/m3] */
#define Vt 565.0                                    /* total liquid volume [m3] = HARDCODED */

#define cell_surf_area 321.0                        /* surface area of cells [m2 / molx] */
#define k_HAcD 3.85e-5                              /* diffusivity of acetic acid [m / h] */
#define pH_EC 5.0                                   /* extracelluar pH  */
#define pKa_HAc 4.756                               /* Acetic acid pKa (source: wikipedia) */
#define Vx 5.89e-5                                  /* bacteria molar volume [m3/molx] */

static double elasticity_matrix_glob[num_rxns][num_compounds];                  /* the elasticity matrix */
static double conc_ref_glob[num_compounds];                                     /* the reference conc [mol/m3] */
static double conc_ini_glob[num_compounds];                                     /* the reference conc [mol/m3] */
static double conc_min_glob[num_compounds];                                     /* the minimum allowable conc [mol/m3] */
static double rates_IC_ref_glob[num_rxns];                                      /* the reference rates [mol/molx/h] */

/* set lagrangian timestep */
#define lagrangian_timestep 1e-4                                                 /* timestep in [s] */
#define eps 1e-15                                                  /* NOT ZERO */

/* for writing */
int writing_timer_part_out = 0;
int writing_timer_part_bin = 0;
int writing_timer_bin_loop = 0;

int timer = 1;

/* for mass transfer */
#define SWITCH_EPS_CORR 1 /* 0 for not taking it into account */



/* FUNCTION to read a certain number of rows and columns from a matrix in a .dat file, specified by filename*/
double *read_matrix(int rows, int cols, const char* filename)
{
    double matrix[rows][cols];
    double *a;
    int i = 0; int j = 0;

    a = ( double * ) malloc ( rows * cols * sizeof ( double ) ); 
    FILE *file_pointer;
    file_pointer = fopen (filename, "r");
    if (file_pointer == NULL)
        return 0;

    for(  i = 0; i < rows; ++i)
    {
        for(  j = 0; j < cols; ++j)
            fscanf(file_pointer, "%lf", &matrix[i][j]);
    }

    for(  i = 0; i < rows; ++i)
    {
        for(  j = 0; j < cols; ++j)
            {
                int offset = i *cols + j;
                a[offset] = matrix[i][j];
            }
    }

    fclose (file_pointer); 
    return a; 
}

/* FUNCTION to read a vector of length n_vals from a .dat file, specified by filename*/
double *read_conc(int num_cmpds, const char* filename)
{
    double conc[num_cmpds];
    double *a;
    int i = 0;

    a = ( double * ) malloc ( num_cmpds * sizeof ( double ) );
    FILE *file_pointer;
    file_pointer = fopen (filename, "r");
    if (file_pointer == NULL)
        return 0;

    for(  i = 0; i < num_cmpds; ++i)
    {
        fscanf(file_pointer, "%lf", &conc[i]);
        a[i] = conc[i];
    }

    fclose (file_pointer); 
    return a; 
}

/* FUNCTION to read a vector of length n_vals from a .dat file, specified by filename*/
double *read_rates(int num_rxn, const char* filename)
{
    double rates[num_rxn];
    double *a;
    int i = 0; 

    a = ( double * ) malloc ( num_rxn * sizeof ( double ) );
    FILE *file_pointer;
    file_pointer = fopen (filename, "r");
    if (file_pointer == NULL)
        return 0;

    for(  i = 0; i < num_rxn; ++i)
    {
        fscanf(file_pointer, "%lf", &rates[i]);
        a[i] = rates[i];    
    }

    fclose (file_pointer); 
    return a; 
}

/* FUNCTION to store the reference rates, concentrations, and elasticity matrix */
DEFINE_ADJUST(get_rates_matrix, domain)
{
    /* only do this when nothing has been read yet (otherwise it is executed too often) */
    if (rates_IC_ref_glob[0] != 0.182)
    {
        /* pointers for concentrations and rates*/
        double *conc_ref_p; 
        double *conc_ini_p; 
        double *conc_min_p; 
        double *rates_IC_ref_p; 
        double *matrix_temp;
        double elasticity_matrix_temp[num_rxns][num_compounds];

        /* read the rates of the reference state [ mol/molx/h ] */
        rates_IC_ref_p = read_rates(num_rxns, "Ref_rates.dat");
        memcpy(rates_IC_ref_glob, rates_IC_ref_p, num_rxns*sizeof(double) ); 


        /* read the concentrations of the reference state [ mol/m3 ]  */
        conc_ref_p = read_conc(num_compounds, "Ref_conc.dat");
        memcpy(conc_ref_glob, conc_ref_p, num_compounds*sizeof(double) ); 


        /* read the initial concentrations [ mol/m3 ]  */
        conc_ini_p = read_conc(num_compounds, "Ini_conc.dat");
        memcpy(conc_ini_glob, conc_ini_p, num_compounds*sizeof(double) ); 


        /* read the minimum allowable concentrations [ mol/m3 ]  */
        conc_min_p = read_conc(num_compounds, "Min_conc.dat");
        memcpy(conc_min_glob, conc_min_p, num_compounds*sizeof(double) ); 


        /*  read the matrix */
        matrix_temp = read_matrix(num_rxns, num_compounds, "eps_matrix.dat");
        for(  i = 0; i < num_rxns; ++i)
        {
            for(  j = 0; j < num_compounds; ++j)
                {
                    int offset = i * num_compounds + j;
                    elasticity_matrix_temp[i][j] = matrix_temp[offset];
                }
        }    
        memcpy(elasticity_matrix_glob, elasticity_matrix_temp, (num_compounds*num_rxns)*sizeof(double) ); 

        
        /* modify the reference concentrations for the case with low ethanol or acetate */
        // double c_Ac_switch = 39.0;
        // if (conc_ini_glob[3] < c_Ac_switch)
        //     {
        //     conc_ref_glob[3] = conc_ini_glob[3];
        //     conc_ref_glob[4] = 0.5;}
        // else
        //     {conc_ref_glob[3] = 0.5;}



        /* clear the memory */
        free(conc_ref_p);
        free(conc_ini_p);
        free(conc_min_p);
        free(rates_IC_ref_p);
        free(matrix_temp);
    }
}

/* FUNCTION to calculate the rates based on the current concentrations - 
no constraints so does not work at super low concentrations -
this is implemented in the conc_derivatives function */
double *rates_fun(double conc_act[num_compounds], double* rates_IC, double rates_IC_ref[num_rxns], double conc_ref[num_compounds], double conc_min[num_compounds], double elasticity_matrix[num_rxns][num_compounds])
{    

    /* definitions */
    double conc_eff[num_compounds];
    double conc_log[num_compounds];
    double rates_ll_cor[num_rxns];
    double r_HAcD = 0; double c_HAc_EC = 0;

    /* check whether the concentrations are too low (c < conc_min)  [ mol/m3 ] */
    for(  i = 0; i < num_compounds; ++i)
    {
        if (conc_act[i] < conc_min[i])
        {conc_eff[i] = conc_min[i];}
        else
        {conc_eff[i] = conc_act[i];}
    }

    /* get the natural logarithm of the IC conc vs its ref [ - ] */
    for(  i = 0; i < num_compounds; ++i)
        conc_log[i] = log(conc_eff[i]/conc_ref[i]);

    /* calculate the linlog corrected rates [ - ] */
    for(  i = 0; i < num_rxns; ++i)
    { 
        double sum_elast = 0;
        for(  j = 0; j < num_compounds; ++j)
        {
            sum_elast += elasticity_matrix[i][j] * conc_log[j];
        }
        
        rates_ll_cor[i] = 1 + sum_elast;
        /* the irreversible reactions (ACAS, AcS, Rnf, Ana, AcX) cannot become negative */
        if (rates_ll_cor[i] < 0 && (i == 1 || i == 4 || i == 8 || i == 9 || i == 10) )
            {rates_ll_cor[i] = 0;}
    }

    /* calculate the total rates of intracellular reactions [ mol/molx/h ]*/
    for(  i = 0; i < num_rxns; ++i)
    { 
        rates_IC[i]= rates_IC_ref[i]*rates_ll_cor[i];                               /* intracellular rates [ mol/mol/h ] */
    }

    // /* Nfn is 3.0 times faster!! */
    // rates_IC[7] = rates_IC[7]*3.0;

    /* calculate the rate of acetate back diffusion [ mol/molx/h ] */
    c_HAc_EC = conc_act[3]*pow(10,pKa_HAc)/(pow(10,pH_EC)+pow(10,pKa_HAc));         /* conc of HAc (extracellular) [ mol/m3l ]*/
    r_HAcD = cell_surf_area*k_HAcD*c_HAc_EC;                                        /* back diffusion rate [ mol/molx/h ] */
    
    rates_IC[11] = r_HAcD;                                                          /* store the back-diffusion as another rate */

    return rates_IC;

}

/* FUNCTION to calculate derivatives */
double *conc_derivatives(double conc_act[num_compounds], double dt, double NP_local, double* dc_dt, double rates_IC_ref[num_rxns], double conc_ref[num_compounds], double conc_min[num_compounds], double elasticity_matrix[num_rxns][num_compounds])
{    

    /* definitions */
    double *rates_IC = ( double * ) calloc( (num_rxns+1) , sizeof(double));
    double conc_eff[num_compounds];                          /* corrected concentrations [mol/m3] */
    double r_HAcD = 0; double mu = 0;                                       /* two rates for easy calculations [mol/molx/h] */

    /* check whether the concentrations are too low (c < conc_min) [ mol/m3 ] */
    for (i = 0; i < num_compounds; ++i){
        if (conc_act[i] < conc_min[i])
            {conc_eff[i] = conc_min[i];}
        else
            {conc_eff[i] = conc_act[i];}
    }

    /* calculate rates [ mol/molx/h ] using the rates function */
    rates_IC = rates_fun(conc_act, rates_IC, rates_IC_ref, conc_ref, conc_min, elasticity_matrix);
    
    /* unpack two of them for easy processing [mol/mol/h] */
    mu = rates_IC[9];                                               /* biomass growth rate */
    r_HAcD = rates_IC[11];                                          /* acetate back-diffusion rate */

        /* dc/dt 0 - 5 are for diffusive compounds:     dcIC/dt = dcEC/dt [mol/m3_L/h] */
    dc_dt[0] = -1.*(rates_IC[1] + rates_IC[2])   * CX ;             /* CO uptake rate [mol/m3/h] (negative for consumption) */
    dc_dt[1] = -2.*(rates_IC[3])                 * CX ;             /* H2 uptake rate [mol/m3/h] (negative for consumption) */
    dc_dt[2] =    (rates_IC[2] - rates_IC[0])   * CX ;              /* CO2 uptake rate [mol/m3/h] (negative for consumption)*/
    dc_dt[3] =    (rates_IC[10] - r_HAcD)       * CX ;              /* Net acetate excretion rate [mol/m3/h] */
    dc_dt[4] =    (rates_IC[5])                 * CX ;              /* EtOH prod rate [mol/m3/h] */
    dc_dt[5] =    (rates_IC[6])                 * CX ;              /* BDO prod rate [mol/m3/h] */
        /* dc/dt 6 - 11 are for compounds that remain in cell [mol/m3_cell/h] */
    dc_dt[6] =  (rates_IC[1] -   rates_IC[4] -    2.*rates_IC[6] - 1./2.*rates_IC[9])               * (1.0/Vx) - mu*conc_act[6];   /* AcCoA balance [mol/m3c/h] */
    dc_dt[7] =  (rates_IC[8] - 2.*rates_IC[1] -     rates_IC[5]  - 1./2.*rates_IC[6] - rates_IC[7]) * (1.0/Vx) - mu*conc_act[7];   /* NADH balance [mol/m3c/h] */
    dc_dt[8] =  (rates_IC[3] + 2.*rates_IC[7] - 1./2.*rates_IC[0]  -     rates_IC[1])               * (1.0/Vx) - mu*conc_act[8];   /* NADPH balance [mol/m3c/h] */
    dc_dt[9] =  (rates_IC[4] +   r_HAcD   -     rates_IC[5]  -     rates_IC[10])                 * (1.0/Vx) - mu*conc_act[9];   /* Ac_i_in balance [mol/m3c/h] */
    dc_dt[10] = (rates_IC[0] -   rates_IC[1])                                                    * (1.0/Vx) - mu*conc_act[10];  /* Formate balance [mol/m3c/h] */
    dc_dt[11] = (-0.5*rates_IC[0] + rates_IC[1] + rates_IC[2] + rates_IC[3] - 
                rates_IC[5] - rates_IC[6] - rates_IC[7] - rates_IC[8])                           * (1.0/Vx) - mu*conc_act[11];  /* Ferredoxin balance [mol/m3c/h] */

    /* --------------- -------------------------------------------------CORRECTION ----------------------------------------------
    
    for the dc/dt when the values became too low            
    check for every compound if the calculated dc/dt term is possible according to the current concentrations and the minimum obtainable ones. 
    In case one or more dc/dt terms are infeasible (dc/dt > dc/dt max), then we have to correct ALL other dc/dt terms based on the one which is the furthest away from the possible dc/dt. 

    We assume that (dcidt/dcjdt)_model = (dcidt/dcjdt)_correction and then dcjdt_correction is the maximum dc/dt of the compound that is furthest away from the possible solution. 
    
    The maximum rate is calculated as the concentration difference between the effective concentration and the minimum concentration, 
    divided by the integration time (dependend on the stage in the RK4 solution). 
    For the extracellular compounds, the total number of particles in that cell is needed (in the 2way algorithm) as the concentration difference has to be shared (are bacteria communists??)
    */
    
    /* definitions */
    double max_dcdt[num_compounds];                                             /* the compound specific max dc/dt (mol/m3/h) */
    double dcdt_ratio[num_compounds];                                           /* stoichiometries */
    double dcdt_diff[num_compounds];                                            /* the difference between the calculated dcdt and its max (should not be negative!) */
    int max_index = 0; double max_rate = 0;                                     /* definitions for maximum determination */
    int limiting_rate = 0;                                                      /* flag to store if some rates becomes negative and should be corrected therefore */
    
    for( i = 0; i < num_compounds; ++i)
    { 
        if (i < 6)
            {max_dcdt[i] = (conc_eff[i]-conc_min[i]*0.9)/(NP_local*dt);}        /* maximum rate EC compounds [mol/m3L/part/h] */
        else
            {max_dcdt[i] = (conc_eff[i]-conc_min[i]*0.9)/dt;}                   /* maximum rate IC compounds [mol/m3c/h] */

        dcdt_diff[i] = dc_dt[i] - max_dcdt[i];                                  /* difference between current rate and its max */

        /* determine the index for the compound with the largest rate difference */
        if (dcdt_diff[i] > max_rate)
            {
                max_rate = dcdt_diff[i];
                max_index = i;
            }
        /* if some dc/dt's cause negative concentrations, assign the flag so that we have to correct the other rates */
        if (dc_dt[i] > max_dcdt[i])                                              
            {limiting_rate = 1;}
    }
    
    /* in case there is one rate bigger than allowed */
    if (limiting_rate == 1)
    {
        /* make the dcdt with the max_index the new reference rate  */
        double dc_dt_ref = max_dcdt[max_index];

        /* correct the other dcdt with that factor, considering the current stoichiometry (dcdt-ratio's)  */
        for(int i = 0; i < num_compounds; ++i)
        {
            dcdt_ratio[i] = dc_dt[i]/(dc_dt[max_index]);          /* original dcdt / dcdt that should be corrected for (use abs value so that the sign remains constant) */
            dc_dt[i] = dc_dt_ref * dcdt_ratio[i];                     /* corrected dcdt */
        }
    }

    free(rates_IC);
    return dc_dt;
}

/* FUNCTION that performs the Runge-Kutta integration */
double *rk4vec ( double t0, double *u0, double dt, double NP_local, double rates_IC_ref[num_rxns], double conc_ref[num_compounds], double conc_min[num_compounds], double elasticity_matrix[num_rxns][num_compounds] )
/*
  Nomenclature

  t = current time
  u = concentration vector (contains concentrations at a time t [mol/m3] )
  dt = lagrangian timestep (in HOURS to be passed to the derivatives function)
  NP_local = number of particles in this grid cell
  f = dc/dt vector

*/
{
    double *f0, *f1, *f2, *f3;
    double *u1, *u2, *u3, *u;
    f0 = (double*) calloc(num_compounds, sizeof(double));
    f1 = (double*) calloc(num_compounds, sizeof(double));
    f2 = (double*) calloc(num_compounds, sizeof(double));
    f3 = (double*) calloc(num_compounds, sizeof(double));
    u1 = (double*) calloc(num_compounds, sizeof(double));
    u2 = (double*) calloc(num_compounds, sizeof(double));
    u3 = (double*) calloc(num_compounds, sizeof(double));
    u = (double*) calloc(num_compounds, sizeof(double));

    /* Get four sample values of the derivative */
    f0 = conc_derivatives (u0, dt, NP_local, f0, rates_IC_ref, conc_ref, conc_min, elasticity_matrix);
    /* 1st update of concentrations */
    for (i = 0; i < num_compounds; ++i )
    {
        u1[i] = u0[i] + dt * f0[i] / 2.0;
    }

    /* 2nd derivative */
    f1 = conc_derivatives (u1, dt/2.0, NP_local, f1, rates_IC_ref, conc_ref, conc_min, elasticity_matrix);
    /* 2nd update of concentrations */
    for (i = 0; i < num_compounds; ++i )
    {
        u2[i] = u0[i] + dt * f1[i] / 2.0;
    }

    /* 3rd derivative */
    f2 = conc_derivatives (u2, dt/2.0, NP_local, f2, rates_IC_ref, conc_ref, conc_min, elasticity_matrix);
    /* 2rd update of concentrations */
    for (i = 0; i < num_compounds; ++i )
    {
        u3[i] = u0[i] + dt * f2[i];
    }

    /* 4th derivative */
    f3 = conc_derivatives (u3, dt, NP_local, f3, rates_IC_ref, conc_ref, conc_min, elasticity_matrix);
    
    /* Combine the derivatives to estimate the solution */
    for (i = 0; i < num_compounds; ++i )
    {
        u[i] = u0[i] + dt * ( f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i] ) / 6.0;

        /* in case the concentration is lower than the minimum allowable one, give the conc that value */
        if ( (u[i] < conc_min[i]) )
            {u[i] = conc_min[i]*0.9;}
    }
    
    free(f0);
    free(f1);
    free(f2);
    free(f3);
    free(u1);
    free(u2);
    free(u3);

    return u;
}



/* FUNCTION: to print the rates + matrix file data. Do this before the run of a simulation to check whether everything went well and all data is loaded */
/* the host_to_node_sync_file command is probably one to the solutions to the parallelization problem */
/*  you can easily execute it from the tui :
/define/user-defined/execute-on-demand "print_rates_matrix::libudf-DPM-2way-80K_MT_Timeav" q */
DEFINE_ON_DEMAND(print_rates_matrix)
{

   /* pointers for concentrations and rates*/
   double *conc_ref_p = calloc(num_compounds, sizeof(double));
   double *conc_min_p = calloc(num_compounds, sizeof(double));
   double *rates_IC_ref_p = calloc(num_compounds, sizeof(double));
   double *matrix_temp = calloc(num_rxns*num_compounds, sizeof(double));
    
    /* storing conc and rates */
    double elasticity_matrix_part[num_rxns][num_compounds];                  /* the elasticity matrix */
    double conc_ref_part[num_compounds];                                     /* the reference conc [mol/m3] */
    double conc_min_part[num_compounds];                                     /* the minimum allowable conc [mol/m3] */
    double rates_IC_ref_part[num_rxns];     

    #if RP_HOST
    int file_ratesICref = host_to_node_sync_file("Ref_rates.dat");
    int file_concref = host_to_node_sync_file("Ref_conc.dat");
    int file_conc_min = host_to_node_sync_file("Min_conc.dat");
    int file_matrix = host_to_node_sync_file("eps_matrix.dat");
    #endif

    #if RP_NODE
    int file_ratesICref = host_to_node_sync_file("Ref_rates.dat");
    int file_concref = host_to_node_sync_file("Ref_conc.dat");
    int file_conc_min = host_to_node_sync_file("Min_conc.dat");
    int file_matrix = host_to_node_sync_file("eps_matrix.dat");
    #endif

    printf("File conc ref = %d bytes at node %d", file_concref, myid);

    /* read the rates of the reference state [ mol/molx/h ] */
    rates_IC_ref_p = read_rates(num_rxns, "Ref_rates.dat");
    for (i = 0; i < num_rxns; ++i)
        rates_IC_ref_part[i] = rates_IC_ref_p[i];
    for (i = 0; i < num_rxns; ++i)
        {
            printf("rates_IC_ref[%d] = %lf \n", i, rates_IC_ref_part[i]);
        }

    /* read the concentrations of the reference state [ mol/m3 ]  */
    conc_ref_p = read_conc(num_compounds, "Ref_conc.dat");
    for (i = 0; i < num_compounds; ++i)
        conc_ref_part[i] = conc_ref_p[i];
    for (i = 0; i < num_compounds; ++i)
    {
        printf("conc_ref_part[%d] = %lf \n", i, conc_ref_part[i]);
    }

    /* read the minimum allowable concentrations [ mol/m3 ]  */
    conc_min_p = read_conc(num_compounds, "Min_conc.dat");
    for (i = 0; i < num_compounds; ++i)
        conc_min_part[i] = conc_min_p[i];
    for (i = 0; i < num_compounds; ++i)
    {
        printf("conc_min_part[%d] = %lf \n", i, conc_min_part[i]);
    }

    /*  read the matrix */
    matrix_temp = read_matrix(num_rxns, num_compounds, "eps_matrix.dat");
    for(  i = 0; i < num_rxns; ++i)
    {
        for(  j = 0; j < num_compounds; ++j)
            {
                int offset = i * num_compounds + j;
                elasticity_matrix_part[i][j] = matrix_temp[offset];
            }
    }    
    /* ... and print the matrix */
    for(  i = 0; i < num_rxns; ++i)
    {
        for(  j = 0; j < num_compounds; ++j)
            {
               printf("matrix[%d][%d] = %lf \n", i,j, elasticity_matrix_part[i][j]);
            }
    }    

    /* clear the memory */
    free(conc_ref_p);
    free(conc_min_p);
    free(rates_IC_ref_p);
    free(matrix_temp);
}


/* FUNCTION to check whether the functions (rate/dcdt/rk4) function properly. 
Does not harm to check it upfront the simulation. Simple using a command in the TUI. 

/define/user-defined/execute-on-demand "print_function_results::libudf-DPM-2way-80K_MT_Timeav" q */
DEFINE_ON_DEMAND(print_function_results)
{
        
    /* first read everything */
    double *conc_ref_p = calloc(num_compounds, sizeof(double));
    double *conc_min_p = calloc(num_compounds, sizeof(double));
    double *rates_IC_ref_p = calloc(num_compounds, sizeof(double));
    double *matrix_temp = calloc(num_rxns*num_compounds, sizeof(double));

    // double *conc_ref_p = aligned_alloc(64, num_compounds * sizeof(double));
    // double *conc_min_p = aligned_alloc(64, num_compounds * sizeof(double));
    // double *rates_IC_ref_p = aligned_alloc(64, num_rxns * sizeof(double));
    // double *matrix_temp = aligned_alloc(64, (num_rxns*num_compounds) * sizeof(double));
    
    /* storing conc and rates */
    double elasticity_matrix_part[num_rxns][num_compounds];                  /* the elasticity matrix */
    double conc_ref_part[num_compounds];                                     /* the reference conc [mol/m3] */
    double conc_min_part[num_compounds];                                     /* the minimum allowable conc [mol/m3] */
    double rates_IC_ref_part[num_rxns];     

    /* read the rates of the reference state [ mol/molx/h ] */
    rates_IC_ref_p = read_rates(num_rxns, "Ref_rates.dat");
    for (i = 0; i < num_rxns; ++i)
        rates_IC_ref_part[i] = rates_IC_ref_p[i];

    /* read the concentrations of the reference state [ mol/m3 ]  */
    conc_ref_p = read_conc(num_compounds, "Ref_conc.dat");
    for (i = 0; i < num_compounds; ++i)
        conc_ref_part[i] = conc_ref_p[i];

    /* read the minimum allowable concentrations [ mol/m3 ]  */
    conc_min_p = read_conc(num_compounds, "Min_conc.dat");
    for (i = 0; i < num_compounds; ++i)
        conc_min_part[i] = conc_min_p[i];

    /*  read the matrix */
    matrix_temp = read_matrix(num_rxns, num_compounds, "eps_matrix.dat");
    for(  i = 0; i < num_rxns; ++i)
    {
        for(  j = 0; j < num_compounds; ++j)
            {
                int offset = i * num_compounds + j;
                elasticity_matrix_part[i][j] = matrix_temp[offset];
            }
    }    

    /* clear the memory */
    free(conc_ref_p);
    free(conc_min_p);
    free(rates_IC_ref_p);
    free(matrix_temp);

    /* now read the data to test Ini_conc is easy but every other set of concentrations can be loaded in an other file 
    - then just adjust the filename. BUT keep in mind while copypaste, otherwise filename is not found >> SIGSEV error */
    double *conc_test = ( double * ) calloc(num_compounds, sizeof(double));
    conc_test = read_conc(num_compounds, "Ini_conc.dat");
    for (i = 0; i < num_compounds; ++i)
        {
            printf("conc_test[%d] = %lf \n", i, conc_test[i]);
        }

    /* get rates */
    double* rates_IC_p = calloc( (num_rxns + 1), sizeof(double));
    rates_IC_p = rates_fun(conc_test, rates_IC_p, rates_IC_ref_part, conc_ref_part, conc_min_part, elasticity_matrix_part);
    for (i = 0; i < num_rxns; ++i)
        {
            printf("rates_IC_[%d] = %lf \n", i, rates_IC_p[i]);
        }

    /* get dcdt */
    double* dcdt = calloc(num_compounds, sizeof(double));
    double dt = (lagrangian_timestep/3600); double NP_local = 1.0;
    dcdt = conc_derivatives (conc_test, dt, NP_local, dcdt, rates_IC_ref_part, conc_ref_part, conc_min_part, elasticity_matrix_part);
    for (i = 0; i < num_compounds; ++i)
        {
            printf("dcdt[%d] = %lf \n", i, dcdt[i]);
        }

    /* get the updated concentrations for 1 timestep */
    double* conc_new = (double*) calloc(num_compounds, sizeof(double));;    /*  Storage vector for the updated conc. */
    conc_new = rk4vec ( 3500.1, conc_test, dt, NP_local, rates_IC_ref_part, conc_ref_part, conc_min_part, elasticity_matrix_part);

    for (i = 0; i < num_compounds; ++i)
        {
            printf("conc_new[%d] = %lf \n", i, conc_new[i]);
        }

    /* free some memory */
    free(conc_test);
    free(rates_IC_p);
    free(dcdt);
    free(conc_new);
  
}


/* FUNCTION to update the particles in every timestep : core of the model */
DEFINE_DPM_SCALAR_UPDATE(part_concentration, c, t, initialize, p)
{    
 
    /* some pointers for cleaning */
    Domain *domain;
    cell_t cell;
    Thread *tr;

    /* cleaner if required: reset the UDMs representing the source terms */
    if (C_UDMI(c, t, 14) > 0)                              /* flag for cleaning the source terms udmis */
    {
        domain = Get_Domain(1);                             /* somehow this has to be here, otherwise it crashes.... */

        /* loop over thread and cells for cleaning them */
        thread_loop_c(tr, domain)
        {
            begin_c_loop(cell, tr)
            {
                C_UDMI(cell, tr, 14) = 0;                   /* flag that we should implement the source terms - new sources have to be calculated so should be 0 */
                C_UDMI(cell, tr, 15) = 0;                   /* two-way flag [1 when we can do two-way coupling] */
                C_UDMI(cell, tr, 18) = 0;                   /* local CO source term from the previous timestep [kg/part/s] */
                C_UDMI(cell, tr, 20) = 0;                   /* local H2 source term from the previous timestep [kg/part/s] */
                C_UDMI(cell, tr, 22) = 0;                   /* local CO2 source term from the previous timestep [kg/part/s] */
                C_UDMI(cell, tr, 86) = 0;                   /* nan-flag */
                C_UDMI(cell, tr, 87) = 0;                   /* stiff-solution-flag */
                C_UDMI(cell, tr, 88) = 0;                   /* check for min conc */
                C_UDMI(cell, tr, 92) = 0;                   /* check for CO uptake */
                C_UDMI(cell, tr, 93) = 0;                   /* check for H2 uptake */
                C_UDMI(cell, tr, 94) = 0;                   /* check for CO uptake */
                C_UDMI(cell, tr, 95) = 0;                   /* check for H2 uptake */

            }
            end_c_loop(cell, tr)
        }
    }

    double TimeStep_part = P_DT(p) ;                                /* in seconds ! */
	int pID = PP_ID(p);

    /* only go through all the stuff for the selected number particles */
    if ( (pID < num_particles) && (C_UDMI(c,t,13) != 0) )
    {

/* ----- LOAD THE CONCENTRATIONS + RATES + MATRIX AGAIN ------------------------------- */
        /* pointers for concentrations and rates*/
        double *conc_ref_p = calloc(num_compounds, sizeof(double));
        double *conc_min_p = calloc(num_compounds, sizeof(double));
        double *rates_IC_ref_p = calloc(num_compounds, sizeof(double));
        double *matrix_temp = calloc(num_rxns*num_compounds, sizeof(double));
        
        /* storing conc and rates */
        double elasticity_matrix_part[num_rxns][num_compounds];                  /* the elasticity matrix */
        double conc_ref_part[num_compounds];                                     /* the reference conc [mol/m3] */
        double conc_min_part[num_compounds];                                     /* the minimum allowable conc [mol/m3] */
        double rates_IC_ref_part[num_rxns];     

        // /* read the rates of the reference state [ mol/molx/h ] */
        rates_IC_ref_p = read_rates(num_rxns, "Ref_rates.dat");
        for (i = 0; i < num_rxns; ++i)
            rates_IC_ref_part[i] = rates_IC_ref_p[i];

        /* read the concentrations of the reference state [ mol/m3 ]  */
        conc_ref_p = read_conc(num_compounds, "Ref_conc.dat");
        for (i = 0; i < num_compounds; ++i)
            conc_ref_part[i] = conc_ref_p[i];

        /* read the minimum allowable concentrations [ mol/m3 ]  */
        conc_min_p = read_conc(num_compounds, "Min_conc.dat");
        for (i = 0; i < num_compounds; ++i)
            conc_min_part[i] = conc_min_p[i];

        /*  read the matrix */
        matrix_temp = read_matrix(num_rxns, num_compounds, "eps_matrix.dat");
        for(  i = 0; i < num_rxns; ++i)
        {
            for(  j = 0; j < num_compounds; ++j)
                {
                    int offset = i * num_compounds + j;
                    elasticity_matrix_part[i][j] = matrix_temp[offset];
                }
        }    

        /* clear the memory */
        free(conc_ref_p);
        free(conc_min_p);
        free(rates_IC_ref_p);
        free(matrix_temp);

    
/* ----- GET SOME DATA FROM THE ENVIRONMENT OF THE PARTICLE -------------------- */

        /* pointer to liquid thread */
        Thread *SUBT_Liq;
        SUBT_Liq = THREAD_SUB_THREAD(t,0);                                         /* set liquid thread */
        
        /* get the number of particles in this cell and its volume */
        double NP_local = MAX(1, C_UDMI(c,t, 16));

        /* get the CO concentration in the liquid phase */
        double clCO, clH2, clCO2;
        clCO = C_UDMI(c,t,1)/MW_CO*1000;                                /* mol/m3 */
        clH2 = C_UDMI(c,t,6)/MW_H2*1000;                                /* mol/m3 */
        clCO2 = C_UDMI(c,t,9)/MW_CO2*1000;                              /* mol/m3 */

        /* allocate particle variables */
        double time, x, y, z;   
        time = PP_TIME(p);                                              /* time of particle [s] (end of timestep...) */                
        x = PP_POS(p)[0];                                               /* x-coordinate [m] (width of reactor) */
        y = PP_POS(p)[1];                                               /* y-coordinate [m] (radial (incl downcomer)) */
        z = PP_POS(p)[2];                                               /* z-coordinate [m] (axial position) */
        
        /* Variables for writing */
        FILE *PartOutput_conc;
        char PATHNAME_Part[250];				                        /* provide folder name */
        char PATHNAME_Part_conc_full[250];			                    /* stores path to output directory particle tracks */
        char cwd[250];                                                  /* storage of current working directory */
        getcwd(cwd, sizeof(cwd));                                       /* string with cwd */
        sprintf(PATHNAME_Part,"%s/Lifelines", cwd);  			        /* folder name */


/* --------- GET THE CURRENT CONCENTRATIONS FROM THE PARTICLE ----------------------  */

        /* Store the concentrations for data processing */
        double conc_act[num_compounds];                                 /* Current concentrations */
        double conc_start[num_compounds];                               /* Conc at start of time step */

        /* define the current concentrations in the particle */
        for(  i = 0; i < num_compounds; ++i)
        {
            conc_act[i] = P_USER_REAL(p, (i+2));                        /* [mol/m3] these are updated during the sub-timestep integration */
            conc_start[i] = P_USER_REAL(p, (i+2));                      /* [mol/m3] stay constant to remember the good old times */
        }
        
        /* adjust the reference concentration for the solventogenesis model 
                The initial values for EtOH does not have to be adjusted since that is an user input in the files */
        
        // double c_Ac_switch = 39.0;
        // if (conc_act[3] < c_Ac_switch)
        //     {
        //     conc_ref_part[3] = conc_act[3];
        //     conc_ref_part[4] = 0.5;}
        // else
        //     {conc_ref_part[3] = 0.5;}

/* --------- SET PARTICLE ENVIRONMENT VARIABLES AS STARTING POINT -------------------- */

        /* ONLY happens at the beginning of an Eulerian timestep */
        if(initialize)
        {
            P_USER_REAL(p,1) = 0;                                       /* reset the counter for writing */
            
            /* Adjust the current concentrations for CO, H2, CO2 */
            conc_act[0] = clCO;                                         /* mol/m3 */
            conc_act[1] = clH2;                                         /* mol/m3 */
            conc_act[2] = clCO2;                                        /* mol/m3 */
        }
        
/* --------  BEGIN OF THE SIMULATION: SET THE INITIAL CONCENTRATIONS -----------===--- */

        /* ONLY happens only at the start of the simulation (particle initialization - flag was needed */
        if(P_USER_REAL(p,0) == 0)                                       /* When the flag is 0, we have to do the initialization */
        {
            /* Read the intracellular concentrations from the data file */
            double* conc_ini = (double*) calloc(num_compounds, sizeof(double)); 
            conc_ini = read_conc(num_compounds, "Ini_conc.dat");
            
            /* But adjust for the local CO, H2, CO2 concentration  */
            conc_ini[0] = clCO;                                         /* mol/m3 */
            conc_ini[1] = clH2;                                         /* mol/m3 */
            conc_ini[2] = clCO2;                                        /* mol/m3 */

            /* initialize the particle scalars */
            P_USER_REAL(p, 0) = 1;                                      /* initialization flag (now not anymore) */
            P_USER_REAL(p, 1) = 0;                                      /* particle writing counter */

            /* overwrite the particle concentrations with the initial conc [mol/m3] 
            P_USER_REAL 2-13 define the particle concentrations 
            set these in the FLUENT DPM MENU */
            for(  i = 0; i < num_compounds; ++i)
            {
                P_USER_REAL(p, (i+2)) = conc_ini[i];                    /* Conc in particle [mol/m3] */
                conc_act[i] = conc_ini[i];                              /* update the actual concentration vector that is used for integation */
            }
            
            /* WRITING if we want so */
            if ((write_per_timestep == 1) && (pID < num_particles_write))
            {
                /* get the dc-dt values here */
                double* dc_dt_start =  (double*) calloc(num_compounds, sizeof(double)); /* aligned_alloc(64, num_compounds * sizeof(double)); */
                dc_dt_start = conc_derivatives(conc_ini, (lagrangian_timestep)/3600.0, NP_local, dc_dt_start, rates_IC_ref_part, conc_ref_part, conc_min_part, elasticity_matrix_part);     /* [mol/m3/h] */

                /* Write data per particle */
                sprintf(PATHNAME_Part_conc_full,"%s/conc_kinmodel_part%d.out",PATHNAME_Part,pID);     /* set the filename for this particular particle */
                PartOutput_conc = fopen(PATHNAME_Part_conc_full,"w");								/* open the file for this particle */
                
                /* write the file for concentration data */
                fprintf(PartOutput_conc,"%d,%f,%f,%f,%f,%f,%f,%f, %f,%f,  %e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e, %e,%e,%e,%e,%e,%e, %e,%e,%e,%e,%e,%e \n", 
                pID, time, x, y, z, clCO, clH2, clCO2,
                P_USER_REAL(p,0), P_USER_REAL(p,1), 
                P_USER_REAL(p,2), P_USER_REAL(p,3), P_USER_REAL(p,4), P_USER_REAL(p,5), P_USER_REAL(p,6), P_USER_REAL(p,7),
                P_USER_REAL(p,8), P_USER_REAL(p,9), P_USER_REAL(p,10), P_USER_REAL(p,11), P_USER_REAL(p,12), P_USER_REAL(p,13), 
                dc_dt_start[0], dc_dt_start[1], dc_dt_start[2], dc_dt_start[3], dc_dt_start[4], dc_dt_start[5], 
                dc_dt_start[6], dc_dt_start[7], dc_dt_start[8], dc_dt_start[9], dc_dt_start[10], dc_dt_start[11]);           

                fclose(PartOutput_conc);

                free(dc_dt_start);
            }
        
            free(conc_ini); 
        }


/* ------ THIS HAPPENS IN EVERY PARTICLE TIMESTEP ------------------------------------------------ 
              update the particle concentrations using sub-particle timestep integration  using the rk4 function */

        /* define time characteristics for integration */
        double int_start_time, int_end_time, max_int_time, current_time;
        int_start_time = PP_TIME(p) - TimeStep_part;                                 /* PP_TIME is the time at the end of the particle timestep, so correct by the timestep duration to go back to start */                                       /* current time of the particle [s] */
        int_end_time = int_start_time + TimeStep_part - lagrangian_timestep;         /* time over which has to be integrated in the while loop (variable per particle!) [s] */
        max_int_time = CURRENT_TIME-lagrangian_timestep;                             /* integrate until the eulerian time has been obtained [s] */
        current_time = int_start_time;                                               /* counter for the while loop [s] */

        /* Now loop over time to do the lagrangian integration for the intracellular compounds. 
        Mind that clCO, clH2, clCO2-extern stay constant, while internally they change */

        /* int_end_time is based on the current particle timestep. Max_int_time is based on Eulerian timestep. It can well be that in the last Particle timestep, the current_time > Eulerian DT 
        and that the particle is then put back in the new timestep at the place it was at the moment the new Eulerian timestep started */
        
        int nan_flag = 0; int stiff_result_flag = 0; int good_result_flag = 0;       /* Some flags for labelling since we like to label things and persons */
        /* let the magic start and loop over time during the particle timestep */
        while( (current_time < (int_end_time)) && (current_time <= max_int_time ))        
        {

            /* calculate concentrations for current timestep */
            double* conc_new = (double*) calloc(num_compounds, sizeof(double));    /*  Storage vector for the updated conc. */
            conc_new = rk4vec ( current_time, conc_act, (lagrangian_timestep/3600), NP_local, rates_IC_ref_part, conc_ref_part, conc_min_part, elasticity_matrix_part);
            
            /* update the time [s] */
            current_time += lagrangian_timestep;

            /* be a bit tolerant and give the solver a new chance to find new concentrations when there are nan values */
            for( i = 0; i < num_compounds; ++i)
            {
                if( fabs(conc_new[i]) != fabs(conc_new[i]) )
                    {
                    conc_new = rk4vec ( current_time, conc_act, (lagrangian_timestep/3600), NP_local, rates_IC_ref_part, conc_ref_part, conc_min_part, elasticity_matrix_part);
                    }
            }

            /* calculate the derivatives one more time for checking the results. The dcdt's are being used to check 
            whether the concentrations are in a good range - to prevent strange drops or so
                        the high and low concentration margin might be 100 times varying from the forward-euler solution */                            
            double* dc_dt_act =  (double*) calloc(num_compounds, sizeof(double)); 
            dc_dt_act = conc_derivatives(conc_act, (lagrangian_timestep)/3600.0, NP_local, dc_dt_act, rates_IC_ref_part, conc_ref_part, conc_min_part, elasticity_matrix_part);     /* [mol/m3/h] */
            double margin_factor = 100.0;                                       
            
            /* now lets do some temporary labelling based on the individual new concetration result */
            int nan_flag_count = 0; int good_result_flag_count = 0; int stiff_result_flag_count = 0; 
            for( i = 0; i < num_compounds; ++i)
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
                for( i = 0; i < num_compounds; ++i)
                {
                    conc_act[i] = conc_start[i];
                    nan_flag = 1;
                }
            }
            else if (good_result_flag_count > 10)
            {
                for( i = 0; i < num_compounds; ++i)
                {
                    conc_act[i] = conc_new[i];
                    good_result_flag = 1;
                }
            }
            else if (stiff_result_flag_count >= 1)
                {
                    stiff_result_flag = 1;
                }
            
            free(dc_dt_act);

            /* particle writing during timestep... */
            if ((write_per_while_timestep == 1) && (pID == 1 ) && (current_time < 3508.5))
            {

                /* get other relevant variables */
                x = PP_POS(p)[0];  
                y = PP_POS(p)[1];
                z = PP_POS(p)[2];

                FILE *PartOutput_conc_234;
                double* dc_dt_start =  (double*) calloc(num_compounds, sizeof(double));
                dc_dt_start = conc_derivatives(conc_act, (lagrangian_timestep)/3600.0, NP_local, dc_dt_start, rates_IC_ref_part, conc_ref_part, conc_min_part, elasticity_matrix_part);     /* [mol/m3/h] */


                /* Write data per particle */
                sprintf(PATHNAME_Part_conc_full,"%s/conc_kinmodel_while_part%d.out",PATHNAME_Part,pID);     /* set the filename for this particular particle */
                PartOutput_conc_234 = fopen(PATHNAME_Part_conc_full,"a");								/* open the file for this particle */

                /* write the file for concentration data */
                fprintf(PartOutput_conc_234," %f,%f,%f,%f,%f,%f,%f,%f,%f,%f, %f,%f,  %e,%e,%e,%e,%e,%e, %e,%e,%e,%e,%e,%e, %d, %e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e \n", 
                current_time, clCO, clH2, clCO2,         conc_new[1], NP_local, rates_IC_ref_part[1],  conc_ref_part[1],  conc_min_part[1], elasticity_matrix_part[0][2],                                                                            /* matlab: 1-8*/
                P_USER_REAL(p,0), P_USER_REAL(p,1),                                                                                         /* matlab: 9-10*/
                conc_act[0],conc_act[1],conc_act[2],conc_act[3],conc_act[4],conc_act[5],
                conc_act[6],conc_act[7],conc_act[8],conc_act[9],conc_act[10],conc_act[11], nan_flag,
                dc_dt_start[0], dc_dt_start[1], dc_dt_start[2], dc_dt_start[3], dc_dt_start[4], dc_dt_start[5],
                dc_dt_start[6], dc_dt_start[7], dc_dt_start[8], dc_dt_start[9], dc_dt_start[10], dc_dt_start[11]);
 
                fclose(PartOutput_conc_234);
                free(dc_dt_start);

                P_USER_REAL(p, 1) += 1;                             /* add the timestep for the writing counter */
            }

            free(conc_new);    
        }


/* ------- FINAL CHECK WHETHER THE LAGRANGIAN INTEGRATION WENT WELL ------------------------------------------*/
/*              with the Two_way_flag we check whether the result is different from the start concentrations 
                it will give permission to do the two-way coupling later on         */      

        /* update the particle data */
        P_USER_REAL(p, 0) = 1;                                          /* initialization number remains 1*/
        
        /* update current concentrations in the particle */
        int conc_min_flag = 0;
        double Two_way_flag = 1.0;                                      /* 1 when we can do the two-way-coupling for this particle */
        for(  i = 0; i < num_compounds; ++i)
        {
            if( isnan(conc_act[i]) == 0  )            /* false when this is NaN - then nothing should happen (something went wrong then...) */
                {
                    P_USER_REAL(p, (i+2)) = conc_act[i];                /* concentrations are stored for next timestep [mol/m3] */                        
                }
            else
                {
                    Two_way_flag = 0; 
                    conc_act[i] = conc_start[i];
                }
            if (conc_act[i] < conc_min_part[i])
            {
                conc_min_flag = 1;
            }                           
        }

/* --------------------------------------- TWO-WAY COUPLING -------------------------------
                                    recalculate + store the q-rates                                 */
        

        /* initialize rate variables to be 0 */
        double dcdt_CO_part = 0; double dcdt_H2_part = 0;
        double rate_part_CO = 0; double rate_part_CO_kg = 0;
        double rate_part_H2 = 0; double rate_part_H2_kg = 0;
        double rate_part_CO2 = 0; double rate_part_CO2_kg = 0;
        
        /* Storage vector for the dcdt terms at start of the particle timestep. 
        The conc differeence in the particle should match this dc/dt for mass balance conservation */
        double* dc_dt_start = (double*) calloc(num_compounds, sizeof(double));
        dc_dt_start = conc_derivatives(conc_start, (lagrangian_timestep)/3600.0, NP_local, dc_dt_start, rates_IC_ref_part, conc_ref_part, conc_min_part, elasticity_matrix_part);     /* [mol/m3/h] */

        /* if the conc. difference is feasilbe, and the integration time is long enough, then calculate the rates for the extracellular space */
        if(Two_way_flag == 1 && ((current_time - int_start_time) > lagrangian_timestep*2) )
        {
            /* some variables for the calculation */
            double timestep_correction = (current_time - int_start_time)/CURRENT_TIMESTEP;          /* since not all Particles have the same timestep, we have to weigh the qCO term by the relative duration */
            double Vol_part = Vt/num_particles;                                                     /* reactor volume per particle [ m3 / part ] */
            double Vol_cell = C_VOLUME(c,t) * C_VOF(c, SUBT_Liq);                                   /* liquid volume in current grid cell [ m3 ] */

            /* calculate the dc/dt-based rate at the timestep start */
            dcdt_CO_part = dc_dt_start[0] / ( 3600.0) * timestep_correction * MW_CO  / 1000.0;      /* [kg/m3/s] */
            dcdt_H2_part = dc_dt_start[1] / ( 3600.0) * timestep_correction * MW_H2  / 1000.0;      /* [kg/m3/s] */
            C_UDMI(c,t, 92) += MAX(-0.1, MIN(0,  dcdt_CO_part * Vol_part / Vol_cell) );                         /* dc-dt based CO uptake: volume integral should be kg/s and approximate the lagrangian uptake */
            C_UDMI(c,t, 93) += MAX(-0.1, MIN(0,  dcdt_H2_part * Vol_part / Vol_cell) );                         /* dc-dt based H2 uptake: volume integral should be kg/s and approximate the lagrangian uptake */

            /* calculate particle-based rates for the extracellular compounds */
            rate_part_CO = (conc_act[0]-conc_start[0])/(current_time - int_start_time);             /* [mol/m3/s]    real CO uptake by particle in the current timestep*/
            rate_part_H2 = (conc_act[1]-conc_start[1])/(current_time - int_start_time);             /* [mol/m3/s] */
            rate_part_CO2 = (conc_act[2]-conc_start[2])/(current_time - int_start_time);            /* [mol/m3/s] */
            rate_part_CO_kg = rate_part_CO * MW_CO  / 1000.0;                                       /* [kg/m3/s]     unit conversion */
            rate_part_H2_kg = rate_part_H2 * MW_H2  / 1000.0;                                       /* [kg/m3/s] */
            rate_part_CO2_kg = rate_part_CO2 * MW_CO2  / 1000.0;                                    /* [kg/m3/s] */
            C_UDMI(c,t, 94) += MIN(0,  rate_part_CO_kg);                                            
            C_UDMI(c,t, 95) += MIN(0,  rate_part_H2_kg);                                            

            /* Assign dc/dt to local cell - used in the source term udf */
            C_UDMI(c,t, 18) += MIN(0, rate_part_CO_kg * Vol_part);                                  /* CO uptake rate (negative) [kg/part/s] */
            C_UDMI(c,t, 20) += MIN(0, rate_part_H2_kg * Vol_part);                                  /* H2 uptake rate (negative) [kg/part/s] */
                /* CO2 rate can either be postive (production) or negative (consumption) */
            if (isnan(rate_part_CO2_kg) == 0)
                {C_UDMI(c,t, 22) += rate_part_CO2_kg * Vol_part;}
            else
                {C_UDMI(c,t, 22) += 0;}

            /* store the two-way flag to check where it went well and where not */
            C_UDMI(c,t, 15) += Two_way_flag;
        }
        else    
        {   
            /* this when we cannot do the two-way coupling */
            C_UDMI(c,t, 18) += 0;
            C_UDMI(c,t, 20) += 0;
            C_UDMI(c,t, 22) += 0;
            Two_way_flag = 0;
        }
        
        /* store if nan or minimum concentration  */
        C_UDMI(c,t, 86) = nan_flag;
        C_UDMI(c,t, 87) = stiff_result_flag;
        C_UDMI(c,t, 88) = conc_min_flag;

/* -------- WRITE THE PARTICLE DATA EACH TIMESTEP ---------------------------------------------- */
        if ((write_per_timestep == 1) && (pID < num_particles_write))
        {
            /* get other relevant variables */
            x = PP_POS(p)[0];  
            y = PP_POS(p)[1];
            z = PP_POS(p)[2];

            /* Write data per particle */
            sprintf(PATHNAME_Part_conc_full,"%s/conc_kinmodel_part%d.out",PATHNAME_Part,pID);     /* set the filename for this particular particle */
            PartOutput_conc = fopen(PATHNAME_Part_conc_full,"a");								/* open the file for this particle */

            /* write the file for concentration data */
            fprintf(PartOutput_conc,"%d,%d, %f,%f, %f,%f,%f,%f,%f,%f, %f,%f,  %e,%e,%e,%e,%e,%e, %e,%e,%e,%e,%e,%e, %f,%d, %e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e, %e,%e,%e,%e \n", 
            pID, myid, TimeStep_part, int_start_time,  x, y, z, clCO, clH2, clCO2,                                                                              /* matlab: 1-8*/
            P_USER_REAL(p,0), P_USER_REAL(p,1),                                                                                         /* matlab: 9-10*/
            P_USER_REAL(p,2), P_USER_REAL(p,3), P_USER_REAL(p,4), P_USER_REAL(p,5), P_USER_REAL(p,6), P_USER_REAL(p,7),                 
            P_USER_REAL(p,8), P_USER_REAL(p,9), P_USER_REAL(p,10), P_USER_REAL(p,11), P_USER_REAL(p,12), P_USER_REAL(p,13),
            Two_way_flag, nan_flag, 
            dc_dt_start[0], dc_dt_start[1], dc_dt_start[2], dc_dt_start[3], dc_dt_start[4], dc_dt_start[5],
            dc_dt_start[6], dc_dt_start[7], dc_dt_start[8], dc_dt_start[9], dc_dt_start[10], dc_dt_start[11], 
            dcdt_CO_part, rate_part_CO_kg, dcdt_H2_part, rate_part_H2_kg);
            fclose(PartOutput_conc);

            P_USER_REAL(p, 1) += 1;                             /* add the timestep for the writing counter */
            
        }

        /* free the memory */
        free(dc_dt_start);             
    }
        /* FINISHED ! */   
    

    /* for the other particles which are not considered in our analysis */
    else if (pID > num_particles)
    {
         /* define the current concentrations in the particle to be 0  */
        for(  i = 0; i < num_compounds; ++i)
        {
            P_USER_REAL(p, (i+2)) = 0;                        /* [mol/m3]  */
        }
    }

}


/* FUNCTION to set the timestep of the particle */
DEFINE_DPM_TIMESTEP(PartTimeStep, p, dt)
{
	double max_time_step = 1e-2;            /* [s] when Euler timestep is smaller, than that one is used */
	if (dt > max_time_step)
	{
		/* change time step to the maximum step size */
		return max_time_step;
	}
	return dt;
}



/* FUNCTION to write lifelines in txt */
DEFINE_EXECUTE_AT_END(write_part_data)
{
    /* write only once per 100 timesteps */
    if (writing_timer_part_out % frequency_write_lifelines == 0)
    {
        /* initializations for loop over injection */
        Injection *I;
        Injection *dpm_injections = Get_dpm_injections();
        Particle *p; 
        int pID;								/* particle ID */

        /* definitions for file writing */
        FILE *PartOutput_conc_tend;
        char PATHNAME_Part[75];				        /* provide folder name */
        char PATHNAME_Part_conc_full[120];			/* stores path to output directory particle tracks */
        char cwd[250];                                                  /* storage of current working directory */
        getcwd(cwd, sizeof(cwd));                                       /* string with cwd */
        sprintf(PATHNAME_Part,"%s/Lifelines", cwd);  			        /* folder name */

        /* and we need pointers for cell and threads  */
        cell_t c;
        Thread *t;	
        
        /* Here we define the variables we want to store  */
        double time, x, y, z;
        double clCO, clH2, clCO2;

        /* Well then, there we go. loop over injections */
        loop(I, dpm_injections)																
        {
            /* Well then, there we go. loop over particles */
            loop(p, I->p)
            {
                /* how do we call our current particle? yes, it's just an integer :)   */
                pID = (int)PP_ID(p);											/* store ID for filename seeking */
            
                if (pID < num_particles_write)
                {

                    /* where is our current particle? In a cell, in a thread */
                    c = PP_CELL(p);
                    t = PP_CELL_THREAD(p);
                    
                    /* load variables of our current particle */
                    time = PP_TIME(p);
                    x = PP_POS(p)[0];  
                    y = PP_POS(p)[1];
                    z = PP_POS(p)[2];
                    clCO = C_UDMI(c,t,1)/MW_CO*1000;                            /* mol/m3 */
                    clH2 = C_UDMI(c,t,6)/MW_H2*1000;                            /* mol/m3 */
                    clCO2 = C_UDMI(c,t,9)/MW_CO2*1000;                         /* mol/m3 */
                    
                    /* get the number of particles in the cell */
                    double NP_local;
                    NP_local = MAX(1, C_UDMI(c,t, 16));

                    double Two_way_flag = C_UDMI(c,t, 15);
                    double nan_flag = C_UDMI(c,t, 86);
                    double stiff_flag = C_UDMI(c,t, 87);

                    /* define the current concentrations in the particle */
                    double conc_act[num_compounds];
                    for(  i = 0; i < num_compounds; ++i)
                    {
                        conc_act[i] = P_USER_REAL(p, (i+2));
                    }
                    /* get the derivatives for checking */
                    double* dc_dt_start = (double*) calloc(num_compounds, sizeof(double));
                    dc_dt_start = conc_derivatives(conc_act, (lagrangian_timestep)/3600.0, NP_local, dc_dt_start, rates_IC_ref_glob, conc_ref_glob, conc_min_glob, elasticity_matrix_glob);     /* [mol/m3/h] */

                    /* Write data per particle */
                    sprintf(PATHNAME_Part_conc_full,"%s/conc_tend_kinmodel_part%d.out",PATHNAME_Part,pID);     /* set the filename for this particular particle */                   
                    /* write the file for concentration data */
                    PartOutput_conc_tend = fopen(PATHNAME_Part_conc_full,"a");								/* open the file for this particle */
                    /* write the file for concentration data */
                    fprintf(PartOutput_conc_tend,"%d,%d, %f,%f,%f,%f,%e,%e,%e, %f,%f,  %e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e, %f, %f,%f,%f, %e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e \n", 
                    pID, myid, time, x, y, z, clCO, clH2, clCO2,                                                                              /* matlab: 1-8*/
                    P_USER_REAL(p,0), P_USER_REAL(p,1),                                                                                         /* matlab: 9-10*/
                    P_USER_REAL(p,2), P_USER_REAL(p,3), P_USER_REAL(p,4), P_USER_REAL(p,5), P_USER_REAL(p,6), P_USER_REAL(p,7),                 
                    P_USER_REAL(p,8), P_USER_REAL(p,9), P_USER_REAL(p,10), P_USER_REAL(p,11), P_USER_REAL(p,12), P_USER_REAL(p,13),
                    NP_local, Two_way_flag, nan_flag, stiff_flag,
                    dc_dt_start[0], dc_dt_start[1], dc_dt_start[2], dc_dt_start[3], dc_dt_start[4], dc_dt_start[5],
                    dc_dt_start[6], dc_dt_start[7], dc_dt_start[8], dc_dt_start[9], dc_dt_start[10], dc_dt_start[11]);
                    fclose(PartOutput_conc_tend);

                    free(dc_dt_start);
                }               
            }
        }
    }
    writing_timer_part_out += 1;
}


/* FUNCTION to write lifelines in binary */
DEFINE_EXECUTE_AT_END(write_part_data_binary)
{
    /* initializations for loop over injection */
	Injection *I;
	Injection *dpm_injections = Get_dpm_injections();
	Particle *p; 
	int pID;								/* particle ID */

	/* definitions for file writing */
    FILE *PartOutput_conc_bin_tend;
    char PATHNAME_Part[250];				/* provide folder name */
    char PATHNAME_Part_conc_full[250];			/* stores path to output directory particle tracks */
    char cwd[250];                                                  /* storage of current working directory */
    getcwd(cwd, sizeof(cwd));                                       /* string with cwd */
    sprintf(PATHNAME_Part,"%s/Lifelines", cwd);  			        /* folder name */

    /* Define struct which will store all data in binary */
    /* make all double as that is easy to read later on */
    typedef struct bin_conc_data_struct{
        double pID;                                     /* Particle ID - col 1 */
        double time, x, y, z, clCO, clH2, clCO2;        /* cols 2 - 8 */
        double part_conc_array[num_compounds];          /* cols 9 - 21 */
        double NP_local;                                    /* col 22 */
        double Two_way_flag, nan_flag, stiff_flag;    /* col 23-24-25 */
        double dcdt_array[num_compounds];               /* cols 26-37 */
    }bin_conc_data_struct;

    /* Initialize global structure */
    bin_conc_data_struct bin_conc_dat;
    
	/* and we need pointers for cell and threads  */
	cell_t c;
	Thread *t;	
	
    /* write only once per 100 timesteps */
    if (writing_timer_part_bin % frequency_write_lifelines == 0)
    {

        /* Well then, there we go. loop over injections */
        loop(I, dpm_injections)																
        {
            /* Well then, there we go. loop over particles */
            loop(p, I->p)
            {
                /* how do we call our current particle?*/
                pID = (int)PP_ID(p);											/* store ID for filename seeking */
            
                if (pID < num_particles_write)
                {

                    /* definitions */
                    bin_conc_dat.pID = (double)PP_ID(p);
                    
                    /* where is our current particle? In a cell, in a thread */
                    c = PP_CELL(p);
                    t = PP_CELL_THREAD(p);
                    
                    /* load variables of our current particle */
                    bin_conc_dat.time = PP_TIME(p);
                    bin_conc_dat.x = PP_POS(p)[0];  
                    bin_conc_dat.y = PP_POS(p)[1];
                    bin_conc_dat.z = PP_POS(p)[2];
                    bin_conc_dat.clCO = C_UDMI(c,t,1)/MW_CO*1000;                            /* mol/m3 */
                    bin_conc_dat.clH2 = C_UDMI(c,t,6)/MW_H2*1000;                            /* mol/m3 */
                    bin_conc_dat.clCO2 = C_UDMI(c,t,9)/MW_CO2*1000;                         /* mol/m3 */

                    /* define the current concentrations in the particle */
                    double conc_act[num_compounds];
                    for(  i = 0; i < num_compounds; ++i)
                    {
                        conc_act[i] = P_USER_REAL(p, (i+2));
                        bin_conc_dat.part_conc_array[i] = P_USER_REAL(p, (i+2));
                    }

                    /* get the number of particles in the cell */
                    double NP_local = MAX(1, C_UDMI(c,t, 16));
                    bin_conc_dat.NP_local = NP_local;
                    
                    bin_conc_dat.Two_way_flag = C_UDMI(c,t, 15);
                    bin_conc_dat.nan_flag = C_UDMI(c,t, 86);
                    bin_conc_dat.stiff_flag = C_UDMI(c,t, 87);
                    
                    /* get the derivatives */
                    double* dc_dt_start = (double*) calloc(num_compounds, sizeof(double));
                    dc_dt_start = conc_derivatives(conc_act, (lagrangian_timestep)/3600.0, NP_local, dc_dt_start, rates_IC_ref_glob, conc_ref_glob, conc_min_glob, elasticity_matrix_glob);     /* [mol/m3/h] */
                    for(  i = 0; i < num_compounds; ++i)
                    {
                        bin_conc_dat.dcdt_array[i] = dc_dt_start[i];
                    }

                    /* Write data per particle */
                    sprintf(PATHNAME_Part_conc_full,"%s/conc_tend_bin_kinmodel_part%d.bin",PATHNAME_Part,pID);     /* set the filename for this particular particle */                

                    /* WRITE IN BINARY */
                    if((PartOutput_conc_bin_tend=fopen(PATHNAME_Part_conc_full,"ab"))==NULL){
                        Message("Error! opening (bin) file\n");
                    }
                    size_t ewb = fwrite(&bin_conc_dat, sizeof(bin_conc_data_struct), 1, PartOutput_conc_bin_tend);
                    if(ewb==0){
                        fclose(PartOutput_conc_bin_tend);
                        Message("Error! writing (bin) file\n");
                    }
                    fclose(PartOutput_conc_bin_tend);

                    free(dc_dt_start);
                }               
            }
        }
    }
    writing_timer_part_bin += 1;
}



/* FUNCTION to write lifelines in binary using a loop in several files per variable
    check overview in which file which variable is since they are numbered */
DEFINE_EXECUTE_AT_END(write_part_big_data_binary_loop)
{

      /* write only once per 100 timesteps */
    if (writing_timer_bin_loop % frequency_write_big == 0)
    {

        /* initializations for loop over injection */
        Injection *I;
        Injection *dpm_injections = Get_dpm_injections();
        Particle *p; 
        int pID;								/* particle ID */

        /* and we need pointers for cell and threads  */
        cell_t c;
        Thread *t;	
        double conc_act[num_compounds];
        int num_variables = num_compounds + 4 + 5;                      /* 12 compounds + 4 location (NP) + 5 q-rates */
        double part_data[num_particles_write_big][ num_variables ];

        /* Well then, there we go. loop over injections */
        loop(I, dpm_injections)																
        {
            /* Well then, there we go. loop over particles */
            loop(p, I->p)
            {
                /* how do we call our current particle?*/
                pID = (int)PP_ID(p);											/* store ID for filename seeking */
            
                if (pID < num_particles_write_big)
                {

                    /* definitions */
                    double *dc_dt_vec = ( double * ) calloc(num_compounds, sizeof(double));                          /* Storage vector for the dcdts. */
                    
                    /* where is our current particle? In a cell, in a thread */
                    c = PP_CELL(p);
                    t = PP_CELL_THREAD(p);

                    /* define the current concentrations in the particle */
                    for(  i = 0; i < num_compounds; ++i)
                    {
                        conc_act[i] = P_USER_REAL(p, (i+2));
                        part_data[pID][i] = conc_act[i];
                    }

                    /* get the number of particles in the cell */
                    double NP_local = MAX(1, C_UDMI(c,t, 16));
                    part_data[pID][num_compounds] = NP_local;

                    /* get the location of the particles */
                    part_data[pID][num_compounds+1] = PP_POS(p)[0];         /* x data */
                    part_data[pID][num_compounds+2] = PP_POS(p)[1];         /* y data */
                    part_data[pID][num_compounds+3] = PP_POS(p)[2];         /* z data */

                    /* calculate the new derivatives */
                    dc_dt_vec = conc_derivatives(conc_act, lagrangian_timestep/3600, NP_local, dc_dt_vec, rates_IC_ref_glob, conc_ref_glob, conc_min_glob, elasticity_matrix_glob);         /* mol/m3 (liq/cell)/h */

                    /* get q-rates */
                    for(  i = 0; i < 5; ++i)
                    {
                        part_data[pID][num_compounds + 4 + i] = dc_dt_vec[i]/CX;        /* qCO - qH2 - qCO2 - qAC - qETOH */
                    }
                    
                    free(dc_dt_vec);
                }               
            }
        }

      
        /* definitions for file writing */
        char PATHNAME_Part[250];				/* provide folder name */
        char cwd[250];                                                  /* storage of current working directory */
        getcwd(cwd, sizeof(cwd));                                       /* string with cwd */
        sprintf(PATHNAME_Part,"%s/Lifelines", cwd);  			        /* folder name */

        char PATHNAME_part_data[250];			                            /* stores path to output directory  */
        int var;                                                            /* index for variables to loop over */
        for( var = 0; var < (num_variables); ++var)
        {
            /* define a struct for variable data storage */         
            typedef struct bin_struct{
            double time;
            double part_data[num_particles_write_big];}bin_struct;
            
            /* Initialize global structure */
            bin_struct bin_data;
            bin_data.time = CURRENT_TIME;

            /* loop over the part_dataset to store the data per particle in the struct */
            for (pID = 0; pID < num_particles_write_big; ++pID )
            {
                bin_data.part_data[pID] = part_data[pID][var];   
            }

            /* prepare for file writing */
            sprintf(PATHNAME_part_data,"%s/Variable_%d.bin",PATHNAME_Part, var);     /* set the filename for the variable */
            FILE *part_data_file;

            if((part_data_file=fopen(PATHNAME_part_data,"ab"))==NULL){
                Message("Error! opening partdata (bin) file\n"); }
            
            int ewb = fwrite(&bin_data, sizeof(bin_struct), 1, part_data_file);
            if(ewb==0){         
                fclose(part_data_file);
                Message("Error! writing part-data (bin) file\n");}
            fclose(part_data_file);

        }
        
    }
    writing_timer_bin_loop += 1;
}


/* FUNCTION to write lifelines in binary using a loop in several files per variable
    check overview in which file which variable is since they are numbered */
DEFINE_EXECUTE_AT_END(write_part_big_data_loop_out)
{

      /* write only once per 100 timesteps */
    if (writing_timer_bin_loop % frequency_write_big == 0)
    {

        /* only do this when nothing has been read yet (otherwise it is executed too often) */
        if (rates_IC_ref_glob[0] != 0.182)
        {
            /* pointers for concentrations and rates*/
            double *conc_ref_p; 
            double *conc_min_p; 
            double *rates_IC_ref_p; 
            double *matrix_temp;
            double elasticity_matrix_temp[num_rxns][num_compounds];

            /* read the rates of the reference state [ mol/molx/h ] */
            rates_IC_ref_p = read_rates(num_rxns, "Ref_rates.dat");
            memcpy(rates_IC_ref_glob, rates_IC_ref_p, num_rxns*sizeof(double) ); 


            /* read the concentrations of the reference state [ mol/m3 ]  */
            conc_ref_p = read_conc(num_compounds, "Ref_conc.dat");
            memcpy(conc_ref_glob, conc_ref_p, num_compounds*sizeof(double) ); 


            /* read the minimum allowable concentrations [ mol/m3 ]  */
            conc_min_p = read_conc(num_compounds, "Min_conc.dat");
            memcpy(conc_min_glob, conc_min_p, num_compounds*sizeof(double) ); 


            /*  read the matrix */
            matrix_temp = read_matrix(num_rxns, num_compounds, "eps_matrix.dat");
            for(  i = 0; i < num_rxns; ++i)
            {
                for(  j = 0; j < num_compounds; ++j)
                    {
                        int offset = i * num_compounds + j;
                        elasticity_matrix_temp[i][j] = matrix_temp[offset];
                    }
            }    
            memcpy(elasticity_matrix_glob, elasticity_matrix_temp, (num_compounds*num_rxns)*sizeof(double) ); 

            /* clear the memory */
            free(conc_ref_p);
            free(conc_min_p);
            free(rates_IC_ref_p);
            free(matrix_temp);
        }


        /* initializations for loop over injection */
        Injection *I;
        Injection *dpm_injections = Get_dpm_injections();
        Particle *p; 
        int pID;								/* particle ID */

        /* and we need pointers for cell and threads  */
        cell_t c;
        Thread *t;	
        int num_variables = num_compounds + 4 + 5;                      /* 12 compounds + 4 location (NP) + 5 q-rates */
        double part_data[num_variables][num_particles_write_big];
        int var;                                                            /* index for variables to loop over */

        for( var = 0; var < (num_variables); ++var)
        {
           
            int pID_ind = 0;
            /* loop over the part_dataset to store initialize as -100 */
            for (pID_ind = 0; pID_ind < num_particles_write_big; ++pID_ind )
            {
               part_data[var][pID_ind] = -100.0;
            }
        }


        /* Well then, there we go. loop over injections */
        loop(I, dpm_injections)																
        {
            /* Well then, there we go. loop over particles */
            loop(p, I->p)
            {
                /* how do we call our current particle?*/
                pID = (int)PP_ID(p);											/* store ID for filename seeking */
            
                if (pID < num_particles_write_big)
                {

                    
                    /* where is our current particle? In a cell, in a thread */
                    c = PP_CELL(p);
                    t = PP_CELL_THREAD(p);

                    double conc_act[num_compounds];
                    /* define the current concentrations in the particle */
                    for(  i = 0; i < num_compounds; ++i)
                    {
                        conc_act[i] = P_USER_REAL(p, (i+2));
                        part_data[i][pID] = conc_act[i];
                    }

                    /* get the number of particles in the cell */
                    double NP_local = MAX(1, C_UDMI(c,t, 16));
                    part_data[num_compounds][pID] = NP_local;

                    /* get the location of the particles */
                    part_data[num_compounds+1][pID] = PP_POS(p)[0];         /* x data */
                    part_data[num_compounds+2][pID] = PP_POS(p)[1];         /* y data */
                    part_data[num_compounds+3][pID] = PP_POS(p)[2];         /* z data */

                    /* calculate the new derivatives */
                    double *dc_dt_vec = ( double * ) calloc(num_compounds, sizeof(double));                          /* Storage vector for the dcdts. */
                    dc_dt_vec = conc_derivatives(conc_act, lagrangian_timestep/3600, NP_local, dc_dt_vec, rates_IC_ref_glob, conc_ref_glob, conc_min_glob, elasticity_matrix_glob);         /* mol/m3 (liq/cell)/h */

                    /* get q-rates */
                    for(  i = 0; i < 5; ++i)
                    {
                        part_data[num_compounds + 4 + i][pID] = dc_dt_vec[i]/CX;        /* qCO - qH2 - qCO2 - qAC - qETOH */
                    }
                    
                    free(dc_dt_vec);
                }               
            }
        }

       
        
        /* definitions for file writing */
        char PATHNAME_Part[250];				/* provide folder name */
        char cwd[250];                                                  /* storage of current working directory */
        getcwd(cwd, sizeof(cwd));                                       /* string with cwd */
        sprintf(PATHNAME_Part,"%s/Lifelines", cwd);  			        /* folder name */

        char PATHNAME_part_data[250];			                            /* stores path to output directory  */      
        double part_data_vec[num_variables][num_particles_write_big];
        double time = CURRENT_TIME;
            
        for( var = 0; var < (num_variables); ++var)
        {
           
            int pID_ind = 0;
            /* loop over the part_dataset to store the data per particle in the struct */
            for (pID_ind = 0; pID_ind < num_particles_write_big; ++pID_ind )
            {
                part_data_vec[var][pID_ind] = PRF_GRHIGH1( part_data[var][pID_ind] );
            }   
        }


        if (I_AM_NODE_ZERO_P)
        {
            for( var = 0; var < (num_variables); ++var)
            {
                int pID_ind = 0;

                /* prepare for file writing */
                sprintf(PATHNAME_part_data,"%s/Variable_%d.out",PATHNAME_Part, var);     /* set the filename for the variable */
                FILE *part_data_file;
                
                part_data_file=fopen(PATHNAME_part_data,"a");
                fprintf(part_data_file, "%f ", time);

                /* loop over the part_dataset to store the data per particle in the struct */
                for (pID_ind = 0; pID_ind < num_particles_write_big; ++pID_ind )
                {
                    fprintf(part_data_file, " %e", part_data_vec[var][pID_ind]) ;
                }
                
                /* end row in part data file */
                fprintf(part_data_file, " \n");

                fclose(part_data_file);

            }
        }
    }
    writing_timer_bin_loop += 1;

}






/* FUNCTION to clean the memory for the number of particles at the end of timestep + put in new results [required] */
DEFINE_EXECUTE_AT_END(UDMI_cleaner)
{
    Domain *domain;
    domain = Get_Domain(1);
    cell_t cell;
    Thread *tr;

/* cleaner if required reset the UDM that counts the particle in each cell and count how many particles can do the the 2way coupling*/

    thread_loop_c(tr, domain)
    {
        begin_c_loop(cell, tr)
        {
            C_UDMI(cell, tr, 16) = 0;              /* local cell = 0 */
        }
        end_c_loop(cell, tr)
    }

    
/* ---------------------- Now update with the current number of particles -------------------- */

    /* initializations for loop over injection */
	Injection *I;
	Injection *dpm_injections = Get_dpm_injections();
	Particle *p; 
	int pID;								/* particle ID */

    /* and we need pointers for cell and threads  */
	cell_t c;
	Thread *t;	

	/* Well then, there we go. loop over injections */
	loop(I, dpm_injections)																
	{
		/* Well then, there we go. loop over particles */
		loop(p, I->p)
		{
			/* how do we call our current particle?*/
			pID = (int)PP_ID(p);											
		
            if (pID < num_particles)
            { 
                /* where is our current particle? In a cell, in a thread */
                c = PP_CELL(p);
                t = PP_CELL_THREAD(p);

                /* update the UDMI */
                C_UDMI(c,t,16) += 1;
            }
        }
    }
}


/* FUNCTION to store some data in the memory [not needed for execution of kinetic model] */
DEFINE_EXECUTE_AT_END(Calculate_source_terms)
{
    /* function to store local particle-based variables to check the influence of the metabolism on the qCO, qEtOH etc */

    /* initializations for loop over injection */
	Injection *I;
	Injection *dpm_injections = Get_dpm_injections();
	Particle *p; 
	int pID;								                                        /* particle ID */

    /* and we need pointers for cell and threads  */
	cell_t c;
	Thread *t;	

    /* only do this when nothing has been read yet (otherwise it is executed too often) */
    if (rates_IC_ref_glob[0] != 0.182)
    {
        /* pointers for concentrations and rates*/
        double *conc_ref_p; 
        double *conc_min_p; 
        double *rates_IC_ref_p; 
        double *matrix_temp;
        double elasticity_matrix_temp[num_rxns][num_compounds];

        /* read the rates of the reference state [ mol/molx/h ] */
        rates_IC_ref_p = read_rates(num_rxns, "Ref_rates.dat");
        memcpy(rates_IC_ref_glob, rates_IC_ref_p, num_rxns*sizeof(double) ); 


        /* read the concentrations of the reference state [ mol/m3 ]  */
        conc_ref_p = read_conc(num_compounds, "Ref_conc.dat");
        memcpy(conc_ref_glob, conc_ref_p, num_compounds*sizeof(double) ); 


        /* read the minimum allowable concentrations [ mol/m3 ]  */
        conc_min_p = read_conc(num_compounds, "Min_conc.dat");
        memcpy(conc_min_glob, conc_min_p, num_compounds*sizeof(double) ); 


        /*  read the matrix */
        matrix_temp = read_matrix(num_rxns, num_compounds, "eps_matrix.dat");
        for(  i = 0; i < num_rxns; ++i)
        {
            for(  j = 0; j < num_compounds; ++j)
                {
                    int offset = i * num_compounds + j;
                    elasticity_matrix_temp[i][j] = matrix_temp[offset];
                }
        }    
        memcpy(elasticity_matrix_glob, elasticity_matrix_temp, (num_compounds*num_rxns)*sizeof(double) ); 

        /* clear the memory */
        free(conc_ref_p);
        free(conc_min_p);
        free(rates_IC_ref_p);
        free(matrix_temp);
    }




    /* definitions */
    double conc_act[num_compounds];

	/* Well then, there we go. loop over injections */
	loop(I, dpm_injections)																
	{
		/* Well then, there we go. loop over particles */
		loop(p, I->p)
		{
			/* how do we call our current particle?*/
			pID = (int)PP_ID(p);											
		
            if (pID < num_particles)
            {
                double *dc_dt_p = calloc(num_compounds, sizeof(double));                           /* Storage vector for the dcdts. */
                double *rates = calloc(num_compounds, sizeof(double));                          /* Storage vector for the rates. */

                /* where is our current particle? In a cell, in a thread */
                c = PP_CELL(p);
                t = PP_CELL_THREAD(p);
                
                /* define the current concentrations in the particle */
                for(  i = 0; i < num_compounds; ++i)
                {
                    conc_act[i] = P_USER_REAL(p, (i+2));
                }

                /* get the number of particles in the cell */
                double NP_local = MAX(1, C_UDMI(c,t, 16));

                /* calculate the new derivatives */
                dc_dt_p = conc_derivatives(conc_act, lagrangian_timestep/3600, NP_local, dc_dt_p, rates_IC_ref_glob, conc_ref_glob, conc_min_glob, elasticity_matrix_glob);         /* mol/m3 (liq/cell)/h */

                /* calculate the current rates */
                rates = rates_fun(conc_act, rates, rates_IC_ref_glob, conc_ref_glob, conc_min_glob, elasticity_matrix_glob);                                                    /* mol/mol/h */

                /* counter the cell for rolling average (i.e. contribution to new average )*/
                C_UDMI(c,t,17) += 1;

                /* calculate for average local concentration / prod rate */
                for(  i = 0; i < num_compounds; ++i)
                {
                    /* new average   =     old average    * contribution to new   + new average/contribution       */
                    C_UDMI(c,t,i+24) = ( C_UDMI(c,t,i+24)*( C_UDMI(c,t,17)-1 ) + conc_act[i] )/C_UDMI(c,t,17);
                }

                /* calculate for the average local update / prod rate */
                for(  i = 0; i < num_compounds; ++i)
                {
                    /* new average   =     old average    * contribution to new   + new average/contribution       */
                    C_UDMI(c,t,i+36) = ( C_UDMI(c,t,i+36)*( C_UDMI(c,t,17)-1 ) + dc_dt_p[i] )/C_UDMI(c,t,17);
                }

                /* calculate the average local IC rates */
                for(  i = 0; i < num_compounds; ++i)
                {
                    /* new average   =     old average    * contribution to new   + new average/contribution       */
                    C_UDMI(c,t,i+48) = ( C_UDMI(c,t,i+48)*( C_UDMI(c,t,17)-1 ) + rates[i] )/C_UDMI(c,t,17);
                }

                /* get q-rates (mol/molx/h) */
                C_UDMI(c, t, 60) = dc_dt_p[0]/CX;       /* CO */
                C_UDMI(c, t, 61) = dc_dt_p[1]/CX;       /* H2 */
                C_UDMI(c, t, 62) = dc_dt_p[2]/CX;       /* CO2 */
                C_UDMI(c, t, 63) = dc_dt_p[3]/CX;       /* Ac- */
                C_UDMI(c, t, 64) = dc_dt_p[4]/CX;       /* EtOH */
                C_UDMI(c, t, 65) = dc_dt_p[5]/CX;       /* BDO */

                /* free memory */
                free(dc_dt_p);
                free(rates);
            }               
        }
    }
}


/* FUNCTION to set the source term for CO in the two way coupling */
DEFINE_SOURCE(CO_Cons_Part, c, t, dS, eqn)
{
    /* definitions */
    double uptake_rate, source;
    double q_CO = C_UDMI(c,t,18);                       /* Arrived is the local CO uptake per particle [kg CO / part / s] */
    double Vc = C_VOLUME(c, t);                         /* Volume of current cell where reaction takes place [m3] */
    double Vl = Vc * C_VOF(c,t);                        /* Liquid volume in cell [m3] */
    double dt = CURRENT_TIMESTEP;                       /* Eulerian timestep [not used] */
    double clCO = C_UDMI(c,t,1);                        /* concentration of CO in the cell [kg/m3] */
    double clCO_min = conc_min_glob[0]*MW_CO/1000;           /* minimum concentration of CO [kg/m3] */

    /* But there is no reaction in headspace! [ so we have to blank that first] */
    double xc[ND_ND]; C_CENTROID(xc, c, t);              /* current locations */
    double Z_coord = xc[2]; double egcor;
    if ((C_VOF(c, t) < 0.3) && (Z_coord > 20))           /* when there is a lot of gas and when we are high in the reactor */
        {egcor = 0;}
    else
        {egcor = 1 - C_VOF(c,t);}                         /* handy to have a correction term for the gas hold-up */
    C_UDMI(c,t,13) = egcor;

    /* calculate the uptake rate based on the local volume and the number of particles */
    uptake_rate = q_CO / Vl;                            /* kg CO / m3 / s */    

    /* no reaction in the headspace */
    if (C_UDMI(c,t,13) == 0)
        {uptake_rate = 0;}                  

    /* if uptake_rate is larger than the amount of dissolved gas present */
    C_UDMI(c,t,90) = 0; 
    if (fabs(uptake_rate) > fabs((clCO - clCO_min)/dt) )
        {uptake_rate = -1*fabs((clCO - clCO_min*0.9)/dt); 
        C_UDMI(c,t,90) = 1; }                               /* flag for CO limitation*/

    /* if uptake_rate rate turns out to be NaN (don't know if possible...) */
    if ( isnan(uptake_rate) == 1 )
        {uptake_rate = 0; 
        C_UDMI(c,t,90) = 1; }                               /* flag for CO limitation*/


    C_UDMI(c,t,19) = uptake_rate;                       /* check whether this all went well */
    C_UDMI(c,t,14) = 1;                                 /* flag for source term cleaning - shall be done during "advacing DPM" */
    

    dS[eqn] = 0.0;
    source = uptake_rate;
    return source;
}


/* FUNCTION to set the source term for H2 in the two way coupling */
DEFINE_SOURCE(H2_Cons_Part, c, t, dS, eqn)
{
    /* definitions */
    double uptake_rate, source;
    double q_H2 = C_UDMI(c,t,20);                       /* Arrived is the local H2 uptake per particle [kg H2 / part / s] */
    double Vc = C_VOLUME(c, t);                         /* Volume of current cell where reaction takes place [m3] */
    double Vl = Vc * C_VOF(c,t);                        /* Liquid volume in cell [m3] */
    double dt = CURRENT_TIMESTEP;
    double clH2 = C_UDMI(c,t,6);                        /* concentration of H2 in the cell [kg/m3] */
    double clH2_min = conc_min_glob[1]*MW_H2/1000;           /* minimum concentration of H2 [kg/m3] */

    /* calculate the uptake rate based on the local volume and the number of particles */
    uptake_rate = q_H2 / Vl ;                          /* kg H2 / m3 /s */
    /* no reaction in headspace! */
    if (C_UDMI(c,t,13) == 0)
        {uptake_rate = 0;}

    /* if uptake_rate is larger than the amount of dissolved gas present */
    C_UDMI(c,t,91) = 0; 
    if (fabs(uptake_rate) > fabs((clH2 - clH2_min)/dt) )
        {uptake_rate = -1*fabs((clH2 - clH2_min*0.9)/dt); 
        C_UDMI(c,t,91) = 1; }                               /* flag for H2 limitation*/

    /* if uptake_rate rate turns out to be NaN (don't know if possible...) */
    if ( isnan(uptake_rate) == 1 )
        {uptake_rate = 0; 
        C_UDMI(c,t,91) = 1; }                               /* flag for H2 limitation*/    
    
    C_UDMI(c,t,21) = uptake_rate;                      /* check whether this all went well */
    
    dS[eqn] = 0.0;
    source = uptake_rate;
    return source;
}

/* FUNCTION to set the source term for CO2 in the two way coupling */
DEFINE_SOURCE(CO2_Cons_Part, c, t, dS, eqn)
{
    /* definitions */
    double uptake_rate, source;
    double q_CO2 = C_UDMI(c,t,22);                      /* Arrived is the local CO2 uptake/prod per particle [kg CO2 / part / s] */
    double Vc = C_VOLUME(c, t);                         /* Volume of current cell where reaction takes place [m3] */
    double Vl = Vc * C_VOF(c,t);                        /* Liquid volume in cell [m3] */
    double dt = CURRENT_TIMESTEP;
    double clCO2 = C_UDMI(c,t,9);                        /* concentration of CO2 in the cell [kg/m3] */
    double clCO2_min = conc_min_glob[2]*MW_CO2/1000;           /* minimum concentration of CO2 [kg/m3] */

    /* calculate the uptake rate based on the local volume and the number of particles */
    uptake_rate = q_CO2 / Vl ;                    /* kg CO2 / m3 /s */
    /* no reaction in headspace! */
    if (C_UDMI(c,t,13) == 0)
        {uptake_rate = 0;}
        
    /* if uptake_rate is larger than the amount of dissolved gas present */
    if ( (uptake_rate < 0)  &&  (fabs(uptake_rate) > fabs((clCO2 - clCO2_min)/dt)) )
        {uptake_rate = -1*fabs((clCO2 - clCO2_min*0.9)/dt); }

    /* if uptake_rate rate turns out to be NaN (don't know if possible...) */
    if ( isnan(uptake_rate) == 1 )
        {uptake_rate = 0; }

    C_UDMI(c,t,23) = uptake_rate;                      /* check whether this all went well */

    dS[eqn] = 0.0;
    source = uptake_rate;
    return source;
}

/* FUNCTION to calculate gradient dependent quantities  */
DEFINE_EXECUTE_AT_END(NP_analysis)
{
    double NP = num_particles;
    double Vp = (double)Vt/NP;                      /* volume per particle [m3] */
    double PI = 3.141592;
   
    Domain *d;					            /* Load domain */
    d = Get_Domain(1);
    cell_t cell;				            /* will loop over all cells c */
    Thread *tr;				            	/* will loop over all threads tr */
    Thread *SUBT_Liq;			            /* sub-thread for liquid */

    /* definitions per iteration */
    double Dt;                              /* turbulent diffusion */
    double Vc;                              /* cell volume */

    /* Definitions for calculations */
    double CCO, R_CO;
    double CH2, R_H2;
    double T_mix_p, T_reac_CO_c, T_reac_H2_c;
    double beta_ms_CO, beta_ms_H2;
    double sigma_p_c;
    double V_term;
    double dt = CURRENT_TIMESTEP;

    /* loop over all cells  */
    thread_loop_c (tr,d)
       {
		
       SUBT_Liq = THREAD_SUB_THREAD(tr,0);                                         /* set liquid thread */

        begin_c_loop_int (cell,tr)
            {

                /* get local flow + cell conditions */
                Dt = C_MU_T(cell,SUBT_Liq)/C_R(cell,SUBT_Liq);                      /* Local turb diffusivity [m2/s] */
                Vc = C_VOLUME(cell,tr);                                             /* local cell volume [m3] */

                /* get local conditions for CO and H2 */
                CCO = C_UDMI(cell, tr, 1);                                          /* local co conc [kg/m3] */
                CH2 = C_UDMI(cell, tr, 6);                                          /* local h2 conc [kg/m3] */
                R_CO = -1*C_UDMI(cell, tr, 19);                                     /* local CO uptake rate due to particle [kg/m3/s] */
                R_H2 = -1*C_UDMI(cell, tr, 21);                                     /* local H2 uptake rate due to particle [kg/m3/s]  */

                if (R_CO < eps)
                   {R_CO = eps;}
                if (R_H2 < eps)
                   {R_H2 = eps;}

                /* calculate reaction based variables */
                T_mix_p =  pow( 3./(4.*PI), (2./3.)) * (1.0/(PI*Dt)) * pow(Vp, (2./3.));        /* mixing time around particle [s] */
                T_reac_CO_c = MIN(1e4, MAX(1e-1, CCO * Vc / (  (R_CO) * Vp)));                  /* CO reaction time in cell [s] - part number being accounted for in rCO */
                T_reac_H2_c = MIN(1e4, MAX(1e-1, CH2 * Vc / (  (R_H2) * Vp)));                  /* H2 reaction time in cell [s] - part number being accounted for in rH2 */
               
                sigma_p_c = pow(NP*Vc/Vt * (1 - Vc/Vt) , (1./2.));                             /* Spatial STD of particle in the cell */

                if(isnan(T_reac_CO_c) == 1)
                    {T_reac_CO_c = 1e4;}
                if(isnan(T_reac_H2_c) == 1)
                    {T_reac_H2_c = 1e4;}

                V_term = 4.*PI*PI/3. * pow( (Vc/Vp), (2./3.));

                /* nothing happens in  headspace! */
                if (C_UDMI(cell,tr,13) == 0)
                    {
                        T_reac_CO_c = 1/eps;
                        T_reac_H2_c = 1/eps;
                        beta_ms_CO = 0;
                        beta_ms_H2 = 0;
                    }
               
                beta_ms_CO = 4.*PI*PI/3. * pow( (Vc/Vp), (2./3.)) * T_mix_p/T_reac_CO_c ;                       /* artificial gradient indicator for CO [-] */
                beta_ms_H2 = 4.*PI*PI/3. * pow((Vc/Vp), (2./3.)) * T_mix_p/T_reac_H2_c;                       /* artificial gradient indicator for H2 [-] */

                if(isnan(beta_ms_CO) == 1)
                    {beta_ms_CO = eps;}
                if(isnan(beta_ms_H2) == 1)
                    {beta_ms_H2 = eps;}

                /* store these local variables in UDMI -> 100 */
                C_UDMI(cell, tr, 81) = T_mix_p;
                C_UDMI(cell, tr, 82) = beta_ms_CO;
                C_UDMI(cell, tr, 83) = beta_ms_H2;                               
                C_UDMI(cell, tr, 84) = sigma_p_c   ;                                                              
                C_UDMI(cell, tr, 85) = V_term;
                C_UDMI(cell, tr, 89) = pow(Vc, (7./6.))/Dt * (1 - pow(Vc/Vt, (1./2.)) );
                /* cell-volume term */

            }
        end_c_loop_int(cell,tr)
        }
}

/* FUNCTION to check whether the data has been read successfully */
DEFINE_ON_DEMAND(UDMI_reset)
{
    Domain *domain;
    domain = Get_Domain(1);
    cell_t cell;
    Thread *tr;
   
   int mem_locs = 93;
/* cleaner if required reset the UDM that counts the particle in each cell */

    thread_loop_c(tr, domain)
    {
        begin_c_loop(cell, tr)
        {
            for (i = 0; i < mem_locs; ++i)
            
            C_UDMI(cell, tr, i) = 0;              /* local cell = 0 */
        }
        end_c_loop(cell, tr)
    }
    Message("UDMI reset done. All UDMI = 0 \n");
}



DEFINE_ADJUST(diss_sum_calculator, d)
{
	double diss_sum_calculated = 0.0;
    Thread *t = Lookup_Thread(d,ZONE_ID);
    Thread *liq = THREAD_SUB_THREAD(t, 0);
    cell_t c;

	/* Get the value of the thread ID from a user-defined Scheme variable */

	/* thread is only used on compute processes */
    
        begin_c_loop_int(c,t)
        {
            diss_sum_calculated += C_D(c,liq)*C_R(c,liq)*C_VOF(c,liq)*C_VOLUME(c,t);
            /* make sure that right thread pointer is used (liquid or mixture...) */
        }
        end_c_loop_int(c,t)
    /*        Message("diss_sum_calculated %f at node %d \n ",diss_sum_calculated, myid); */

    diss_sum_calculated = PRF_GRSUM1(diss_sum_calculated);


    begin_c_loop_int (c, t)
    {
        C_UDMI(c, t, 12) = diss_sum_calculated;
    }
    end_c_loop_int (c, t)

}





/* CO mass transfer */
DEFINE_MASS_TRANSFER(CO_gas_liquid_MT,cell,thread,from_index,from_species_index, to_index, to_species_index)

{
/* Define gas & liquid sub-threads*/
	Thread *gas = THREAD_SUB_THREAD(thread, from_index);
	Thread *liq = THREAD_SUB_THREAD(thread, to_index);

/* Define Variables */
	real m_lg_CO;
	real y_CO;
    real c_CO;
	real c_s_CO;
	real a;
    real k_L_CO;
	real k_L_CO_LS;
    real k_L_CO_higbie;
    double f_cor;
    double diss_sum;

/*  Local conditions */
    real eG =  C_VOF(cell,gas);         /* m3G/m3R, local gas hold-up */
	real P = C_P(cell,thread) + 101325; /* Pa, Reactor Pressure*/
	real dissipation = C_D(cell,liq);		    /* Turbulence energy dissipation, m2/s3 */
	real w_CO = C_YI(cell,gas,0);       /* CO mass fraction */
	real w_H2 = C_YI(cell, gas, 1);     /* H2 mass fraction */
	real w_CO2 = C_YI(cell,gas,2);      /* CO2 mass fraction */

/* Liquid properties */
    double density_L = C_R(cell,liq);     /* kg/m3, Liquid Density */
    real mu_l = C_MU_L(cell,liq);	    /* Pa.s */
    real nu_l = mu_l/density_L;		    /* m2/s */

/* Gas properties */
    real d_b = C_PHASE_DIAMETER(cell,gas);  /* Bubble diameter, m */

/* Diffusion coefficients (see LST master thesis) */
	double D_L_CO = 2.71*pow(10,-9); /* m^2/s, Diffusion coefficient at 37C*/
	
/* Henry's coefficient: (see excel sheet for unit conversion) at 37C  */
	double H_CO = 2.30*pow(10,-7); /* kg/m^3/Pa */

/* Calculating a [interfacial area per unit volume] (Assuming spherical bubbles) */
	a = 6 * eG / d_b;                   /* 1/m */

/* Calculating MT Coefficient using corrected LS relation*/
    double Pmin = 506.7;                  /* Minimum power input W/m3, Roels Heijnen 1980 */
    double Ptot = Pmin * 565;           /* Total power input W */

    if(SWITCH_EPS_CORR == 0)
        {f_cor = 1.0;}
    else
        {
            diss_sum = C_UDMI(cell,thread,12); 
            f_cor = Ptot/diss_sum;          
        }

    k_L_CO_LS =  0.4 * pow(D_L_CO,0.5) * pow((dissipation*f_cor) / nu_l,0.25); /* m/s */
/* Calculating MT Coefficient using Higbie relation*/
	real u_slip = C_U(cell,gas) - C_U(cell,liq); /* m/s */
	real v_slip = C_V(cell,gas) - C_V(cell,liq); /* m/s */
	real w_slip = C_W(cell,gas) - C_W(cell,liq); /* m/s */
	real vel_slip = pow(pow(u_slip,2) + pow(v_slip,2) + pow(w_slip,2),0.5); /* m/s */
    k_L_CO_higbie = pow(4 * D_L_CO * vel_slip / d_b / 3.14,0.5); /* m/s */

/* Pick maximum kL value */
    if (k_L_CO_LS > k_L_CO_higbie)
        {k_L_CO = k_L_CO_LS;}
        else
        {k_L_CO = k_L_CO_higbie;}

/* Limiting k_L [(1) if volume fraction > 0.4] [(2) if k_L > 0.01] */
	if (C_VOF(cell,gas) >= 0.3)
		{k_L_CO = 0;}
	if (k_L_CO > 0.01)
		{k_L_CO = 0.01;}

/* Calculating concentrations */
	y_CO = (w_CO/MW_CO)/((w_CO/MW_CO)+(w_CO2/MW_CO2)+(w_H2/MW_H2));  /* -, CO gas mmole fraction*/
	c_s_CO = H_CO * P * y_CO;               /* kg/m^3 */
	c_CO = density_L * C_YI(cell,liq,0);    /* kg/m^3 */

/* Calculating Mass Transfer Rate */
	m_lg_CO = k_L_CO * a * (c_s_CO - c_CO); /* kg/m^3/s */

/* Storing check values */
	C_UDMI(cell,thread,0) = c_s_CO;
	C_UDMI(cell,thread,1) = c_CO;
	C_UDMI(cell,thread,2) = m_lg_CO;
	C_UDMI(cell,thread,3) = k_L_CO * a;
	C_UDMI(cell,thread,4) = k_L_CO ;
return (m_lg_CO);
}



/* H2 mass transfer */
DEFINE_MASS_TRANSFER(H2_gas_liquid_MT,cell,thread,from_index,from_species_index, to_index, to_species_index)

{
/* Define gas & liquid sub-threads*/
	Thread *gas = THREAD_SUB_THREAD(thread, from_index);
	Thread *liq = THREAD_SUB_THREAD(thread, to_index);

/* Define Variables */
	real m_lg_H2;
	real y_H2;
    real c_H2;
	real c_s_H2;
	real a;
    real k_L_H2;
	real k_L_H2_LS;
    real k_L_H2_higbie;	
    double f_cor;
    double diss_sum;

/*  Local conditions */
    real eG =  C_VOF(cell,gas);         /* m3G/m3R, local gas hold-up */
	real P = C_P(cell,thread) + 101325; /* Pa, Reactor Pressure*/
	real dissipation = C_D(cell,liq);		    /* Turbulence energy dissipation, m2/s3 */
	real w_CO = C_YI(cell,gas,0);       /* CO mass fraction */
	real w_H2 = C_YI(cell, gas, 1);     /* H2 mass fraction */
	real w_CO2 = C_YI(cell,gas,2);      /* CO2 mass fraction */

/* Liquid properties */
    double density_L = C_R(cell,liq);     /* kg/m3, Liquid Density */
    real mu_l = C_MU_L(cell,liq);	    /* Pa.s */
    real nu_l = mu_l/density_L;		    /* m2/s */
/* Gas properties */
    real d_b = C_PHASE_DIAMETER(cell,gas);  /* Bubble diameter, m */
/* Diffusion coefficients (see LST master thesis) */
	double D_L_H2 = 6.01*pow(10,-9); /* m^2/s, Diffusion coefficient at 37 C*/

/* Henry's coefficient: (see excel sheet for unit conversion) at 37C  */
	double H_H2 = 1.47*pow(10,-8); /* kg/m^3/Pa */

/* Calculating a [interfacial area per unit volume] (Assuming spherical bubbles) */
	a = 6 * eG / d_b;                   /* 1/m */

/* Calculating MT Coefficient */
    double Pmin = 506.7;                  /* Minimum power input W/m3, Roels Heijnen 1980 */
    double Ptot = Pmin * 565;           /* Total power input W */

    if(SWITCH_EPS_CORR == 0)
        {f_cor = 1.0;}
    else
        {
            diss_sum = C_UDMI(cell,thread,12); 
            f_cor = Ptot/diss_sum;          
        }
	k_L_H2_LS =  0.4 * pow(D_L_H2,0.5) * pow((dissipation*f_cor) / nu_l,0.25); /* m/s */

/* Calculating MT Coefficient using Higbie relation*/
	real u_slip = C_U(cell,gas) - C_U(cell,liq); /* m/s */
	real v_slip = C_V(cell,gas) - C_V(cell,liq); /* m/s */
	real w_slip = C_W(cell,gas) - C_W(cell,liq); /* m/s */
	real vel_slip = pow(pow(u_slip,2) + pow(v_slip,2) + pow(w_slip,2),0.5); /* m/s */
    k_L_H2_higbie = pow(4 * D_L_H2 * vel_slip / d_b / 3.14,0.5); /* m/s */

/* Pick maximum kL value */
    if (k_L_H2_LS > k_L_H2_higbie)
        {k_L_H2 = k_L_H2_LS;}
        else
        {k_L_H2 = k_L_H2_higbie;}

/* Limiting k_L [(1) if volume fraction > 0.4] [(2) if k_L > 0.01] */
	if (C_VOF(cell,gas) >= 0.3)
		{k_L_H2 = 0;}
	if (k_L_H2 > 0.01)
		{k_L_H2 = 0.01;}

/* Calculating concentrations */
	y_H2 = (w_H2/MW_H2)/((w_CO/MW_CO)+(w_CO2/MW_CO2)+(w_H2/MW_H2));  /* -, CO gas mmole fraction*/
	c_s_H2 = H_H2 * P * y_H2;               /* kg/m^3 */
	c_H2 = density_L * C_YI(cell,liq,1);    /* kg/m^3 */
/* Calculating Mass Transfer Rate */
	m_lg_H2 = k_L_H2 * a * (c_s_H2 - c_H2); /* kg/m^3/s */

/* Storing check values */
	C_UDMI(cell,thread,5) = c_s_H2;
	C_UDMI(cell,thread,6) = c_H2;
	C_UDMI(cell,thread,7) = m_lg_H2;
return (m_lg_H2);
}



/* CO2 mass transfer */
DEFINE_MASS_TRANSFER(CO2_gas_liquid_MT,cell,thread,from_index,from_species_index, to_index, to_species_index)

{
/* Define gas & liquid sub-threads*/
	Thread *gas = THREAD_SUB_THREAD(thread, from_index);
	Thread *liq = THREAD_SUB_THREAD(thread, to_index);

/* Define Variables */
	real m_lg_CO2;
	real y_CO2;
    real c_CO2;
	real c_s_CO2;
	real a;
    real k_L_CO2;
	real k_L_CO2_LS;
    real k_L_CO2_higbie;
    double f_cor;	
    double diss_sum;

/*  Local conditions */
    real eG =  C_VOF(cell,gas);         /* m3G/m3R, local gas hold-up */
	real P = C_P(cell,thread) + 101325; /* Pa, Reactor Pressure*/
	real dissipation = C_D(cell,liq);		    /* Turbulence energy dissipation, m2/s3 */
	real w_CO = C_YI(cell,gas,0);       /* CO mass fraction */
	real w_H2 = C_YI(cell, gas, 1);     /* H2 mass fraction */
	real w_CO2 = C_YI(cell,gas,2);      /* CO2 mass fraction */

/* Liquid properties */
    double density_L = C_R(cell,liq);     /* kg/m3, Liquid Density */
    real mu_l = C_MU_L(cell,liq);	    /* Pa.s */
    real nu_l = mu_l/density_L;		    /* m2/s */
/* Gas properties */
    real d_b = C_PHASE_DIAMETER(cell,gas);  /* Bubble diameter, m */
/* Diffusion coefficients (see LST master thesis) */
	double D_L_CO2 = 2.56*pow(10,-9); /* m^2/s, Diffusion coefficient of CO2 at 37C */

	/* Henry's coefficient: (see excel sheet for unit conversion) at 37C  */
	double H_CO2 = 1.06*pow(10,-5); /* kg/m^3/Pa at 37C */

/* Calculating a [interfacial area per unit volume] (Assuming spherical bubbles) */
	a = 6 * eG / d_b;                   /* 1/m */

/* Calculating MT Coefficient */
    double Pmin = 506.7;                  /* Minimum power input W/m3, Roels Heijnen 1980 */
    double Ptot = Pmin * 565;           /* Total power input W */

    if(SWITCH_EPS_CORR == 0)
        {f_cor = 1.0;}
    else
        {
            diss_sum = C_UDMI(cell,thread,12); 
            f_cor = Ptot/diss_sum;          
        }	
        k_L_CO2_LS =  0.4 * pow(D_L_CO2,0.5) * pow((dissipation*f_cor) / nu_l,0.25); /* m/s */

/* Calculating MT Coefficient using Higbie relation*/
	real u_slip = C_U(cell,gas) - C_U(cell,liq); /* m/s */
	real v_slip = C_V(cell,gas) - C_V(cell,liq); /* m/s */
	real w_slip = C_W(cell,gas) - C_W(cell,liq); /* m/s */
	real vel_slip = pow(pow(u_slip,2) + pow(v_slip,2) + pow(w_slip,2),0.5); /* m/s */
    k_L_CO2_higbie = pow(4 * D_L_CO2 * vel_slip / d_b / 3.14,0.5); /* m/s */

/* Pick maximum kL value */
    if (k_L_CO2_LS > k_L_CO2_higbie)
        {k_L_CO2 = k_L_CO2_LS;}            
        else
        {k_L_CO2 = k_L_CO2_higbie;}             

/* Limiting k_L [(1) if volume fraction > 0.4] [(2) if k_L > 0.01] */
	if (C_VOF(cell,gas) >= 0.3)
		{k_L_CO2 = 0;}
	if (k_L_CO2 > 0.01)
		{k_L_CO2 = 0.01;}

/* Calculating concentrations */
	y_CO2 = (w_CO2/MW_CO2)/((w_CO/MW_CO)+(w_CO2/MW_CO2)+(w_H2/MW_H2));  /* -, CO gas mmole fraction*/
	c_s_CO2 = H_CO2 * P * y_CO2;               /* kg/m^3 */
	c_CO2 = density_L * C_YI(cell,liq,2);    /* kg/m^3 */
/* Calculating Mass Transfer Rate */
	m_lg_CO2 = k_L_CO2 * a * (c_s_CO2 - c_CO2); /* kg/m^3/s */

/* Storing check values */
	C_UDMI(cell,thread,8) = c_s_CO2;
	C_UDMI(cell,thread,9) = c_CO2;
	C_UDMI(cell,thread,10) = m_lg_CO2;
    C_UDMI(cell,thread,11) = f_cor;

return (m_lg_CO2);
}


/* executes at end of every timestep */
DEFINE_EXECUTE_AT_END(Time_average_std)
{
  Domain *d;					/* Load domain */
  cell_t cell;				/* will loop over all cells c */
  Thread *tr;					/* will loop over all threads tr */
  double BU, M2;				/* backup variable for computation */
  real cCO, cH2, cCO2;
  real qCO, qH2, qCO2, qETOH;

  /* if multiphase specify liquid and gas thread. */ 
  Thread *SUBT_Gas;			/* sub-thread for gas */
  d = Get_Domain(1);

  /* loop over all cells to get average and std data per cell */
  thread_loop_c (tr,d)
  {

    SUBT_Gas = THREAD_SUB_THREAD(tr,1); /* set gas thread */

    /* this list can be made arbitrarily long. Make sure to refer to right thread (for epsilon, 
    in mixture k-epsilon mode this is the mixture level thread tr, in dispersed mode this is the liquid thread)
    make sure to include UDM in simulation! */

    begin_c_loop (cell,tr)
    {

      BU = C_UDMI(cell,tr,66);
      C_UDMI(cell,tr,66) = (C_UDMI(cell,tr,66)*(timer-1) + C_VOF(cell,SUBT_Gas))/(timer) ;    		      /* Gas fraction */

      /* conc_averages */
      BU = C_UDMI(cell,tr,67);
      cCO = C_UDMI(cell,tr,1);
      C_UDMI(cell,tr,67) = (C_UDMI(cell,tr,67)*(timer-1) + cCO )/(timer) ;    			          /* CO conc (g/L) - average */
      M2 = C_UDMI(cell,tr,68)+(cCO - BU)*(cCO - C_UDMI(cell,tr,67));	                    		/* CO conc (g/L) - variance */
      C_UDMI(cell,tr,68) = M2/timer;

      BU = C_UDMI(cell,tr,69);
      cH2 = C_UDMI(cell,tr, 6);
      C_UDMI(cell,tr,69) = (C_UDMI(cell,tr,69)*(timer-1) + cH2)/(timer) ;    			            /* H2 conc (g/L) - average */
      M2 = C_UDMI(cell,tr,70)+(cH2 - BU)*(cH2 - C_UDMI(cell,tr,69)); 	                        /* H2 conc (g/L) - variance */
      C_UDMI(cell,tr,70) = M2/timer;

      BU = C_UDMI(cell,tr,71);
      cCO2 = C_UDMI(cell,tr, 9);
      C_UDMI(cell,tr,71) = (C_UDMI(cell,tr,71)*(timer-1) + cCO2 )/(timer) ;    		            /* CO2 conc (g/L) - average */
      M2 = C_UDMI(cell,tr,72)+(( cCO2 - BU)*( cCO2 - C_UDMI(cell,tr,71)));                    /* CO2 conc (g/L) - variance */
      C_UDMI(cell,tr,72) = M2/timer;

      /* q-averages */
      BU = C_UDMI(cell,tr,73);
      qCO = C_UDMI(cell,tr, 60);
      C_UDMI(cell,tr,73) = (C_UDMI(cell,tr,73)*(timer-1) + qCO )/(timer) ;    		            /* CO q rate (mol/mol/h) - average */

      BU = C_UDMI(cell,tr,74);
      qH2 = C_UDMI(cell,tr, 61);
      C_UDMI(cell,tr,74) = (C_UDMI(cell,tr,74)*(timer-1) + qH2)/(timer) ;    			            /* H2 q rate (mol/mol/h)  - average */

      BU = C_UDMI(cell,tr,75);
      qCO2 = C_UDMI(cell,tr, 62);
      C_UDMI(cell,tr,75) = (C_UDMI(cell,tr,75)*(timer-1) + qCO2 )/(timer) ;    		            /* CO2q rate (mol/mol/h)  - average */

      BU = C_UDMI(cell,tr,76);
      qETOH = C_UDMI(cell,tr, 64);
      C_UDMI(cell,tr,76) = (C_UDMI(cell,tr,76)*(timer-1) + qETOH )/(timer) ;    		            /* ETOH rate (mol/mol/h)  - average */
      M2 = C_UDMI(cell,tr,77)+(qETOH - BU)*(qETOH - C_UDMI(cell,tr,76));	                    		/* qEtOH - variance */
      C_UDMI(cell,tr,77) = M2/timer;

      /* COV FOR CO AND H2 */
      if (C_UDMI(cell,tr,13) != 0)
      {
        if( (C_UDMI(cell, tr, 68) != 0) && (C_UDMI(cell, tr, 67) != 0) )
            {C_UDMI(cell,tr, 78) = ( pow(C_UDMI(cell,tr, 68), (1./2.))  )/C_UDMI(cell,tr, 67) ;}      /* COV FOR CL-CO */

        if( (C_UDMI(cell, tr, 70) != 0) && (C_UDMI(cell, tr, 69) != 0) )
          {C_UDMI(cell,tr, 79) = ( pow(C_UDMI(cell,tr, 70), (1./2.))  )/C_UDMI(cell,tr, 69) ;}        /* COV FOR CL-H2 */

        if( (C_UDMI(cell, tr, 77) != 0) && (C_UDMI(cell, tr, 76) != 0) )
          {C_UDMI(cell,tr, 80) = ( pow(C_UDMI(cell,tr, 77), (1./2.))  )/fabs(C_UDMI(cell,tr, 76)) ;}      /* COV FOR Q-ETOH */
      }
      else
      {
        C_UDMI(cell,tr, 78) = 1e-15 ;      /* COV FOR CL-CO */
        C_UDMI(cell,tr, 79) = 1e-15 ;      /* COV FOR CL-H2 */
        C_UDMI(cell,tr, 80) = 1e-15 ;
      }
    if(C_UDMI(cell,tr,76) == 0)
        {C_UDMI(cell,tr, 80) = 1;}

    }	
    end_c_loop (cell,tr)
  
  }
  
  timer = timer+1;	
}