/*function to compute carbon monoxide uptake (CO)*/

# include "udf.h"
# include "mem.h"
# include "sg_mphase.h"


/* Define global constants */
real MW_CO = 28.01;		    /* Molecular mass: g/mol */
real MW_CO2 = 44.01;		/* Molecular mass: g/mol */
real MW_H2 = 2.016;		    /* Molecular mass: g/mol */
real MW_x = 24.6;			/* Molecular (bio)mass: g/mol */
real c_X = 25; 			    /* Biomass concentration: kg_X/m3*/

/* Global variables */
real q_CO;			        /* Spec. Uptake Rate: kg_CO/kgx/s */
real r_CO;			        /* Uptake Rate: kg_CO/m3/s */
real q_H2;                  /* Uptake Rate: kg_H2/m3/s */
real r_H2;                  /* Uptake Rate: kg_H2/m3/s */

/* executes at end of every time step*/
DEFINE_SOURCE(CO_UPTAKE, cell, thread, dS, eqn)
{
    /* Define Variables */
    real source;		                /* Source Term: kg_CO/m3/s */
    real c_CO; 			                /* CO concentration (liquid): kg/m3 */
    real q_CO_max; 		                /* Maximum uptake rate: kg_CO/kg_X/s */
    real K_M, K_I; 		                /* Monod constants: kg_CO/m3 */
    real density_L = C_R(cell,thread);  /* Density: kg/m3 */
 
    /* Calculating Variables */

    /* Maximum uptake rate: */
    q_CO_max = 1.459/3600;				/* molCO/molx/s */
    q_CO_max = q_CO_max / MW_x;			/* mol CO / gx /h*/

    /* Monod Constants: mM -> kg_CO/m3 */
    K_M = 0.042*MW_CO/1000;		        /* Monod satutration const: mM -> kg_CO/m3 */
    K_I = 0.246 * MW_CO/1000;		    /* Monod Inhibition const: mM -> kg_CO/m3 */

    /* Calculating concentration of CO */
    /* CO concentration (liquid): kg/m3 */
    c_CO = density_L * C_YI(cell,thread,0);

    /* Calculating uptake Rate */
    /* Specific uptake rate: kmol_CO/kg_X/s */
    q_CO = (q_CO_max * c_CO)/ (c_CO + K_M + (pow(c_CO,2)/K_I));

    /* Limiting Biomass Uptake in Droplets and Spray */
    if (C_VOF(cell,thread) <= 0.25)
        {q_CO = 0;}

    /* Defining CO uptake Rate */
    /* Uptake Rate: kmol_CO/m3/s */
    r_CO = c_X * q_CO;

    /* Output source terms: kg_CO/m3/s */
    source = -r_CO * MW_CO;

    /* Storing check values */
    C_UDMI(cell,thread,22) = q_CO*MW_x*3600;			/* Specific uptake rate: kmol_CO/kg_X/s */
    C_UDMI(cell,thread,23) = source;		            /* Uptake Rate: kg_CO/m3/s */ 

    /* Returning Source Term for */	
    dS[eqn]=0.0;
    return source;
}




/* executes at end of every time step*/
DEFINE_SOURCE(H2_UPTAKE, cell, thread, dS, eqn)
{
    /* Define Variables */
    real source;		                    /* Source Term: kg_H2/m3/s */
    real c_H2; 		                        /* H2 concentration (liquid): kg/m3 */
    real c_CO; 		                        /* CO concentration (liquid): kg/m3 */
    real q_H2_max; 		                    /* Maximum uptake rate: kg_H2/kg_X/s */
    real K_M; 		                        /* Monod constants: kg_H2/m3 */
    real density_L = C_R(cell,thread);      /* Density: kg/m3 */
    real K_I_CO, I_CO;                      /* CO inhibiton */

    /* Calculating Variables */
    /* Maximum uptake rate: */
    q_H2_max = 2.565/3600;					/* mol H2/molx/s*/
    q_H2_max = q_H2_max / MW_x;			    /* mol H2 / gx /h*/

    /* Monod Constants: mM -> kg_CO/m3 */
    K_M = 0.025*MW_H2/1000;		            /* Monod satutration const: mM -> kg_H2/m3 */

    /* Calculating concentration of CO and H2 */
    c_CO = density_L * C_YI(cell,thread,0); /* Concentration in liquid (kg/m3) */
    c_H2 = density_L * C_YI(cell,thread,1); /* Concentration in liquid (kg/m3) */

    /* Calculating inhibition by CO */
    /* Acetate concentration (liquid): kg/m3 */
    K_I_CO = 0.025 * MW_CO / 1000;			/* Inhibition constant (kg/m3) */
    I_CO = 1 / (1 + c_CO / K_I_CO);			/* Inhibition term (-) */

    /* Calculating uptake Rate */
    q_H2 = (q_H2_max * c_H2)/ (c_H2 + K_M) * I_CO;      /* Specific uptake rate: kmol_H2/kg_X/s */

    /* Limiting Biomass Uptake in Droplets and Spray */
    if (C_VOF(cell,thread) <= 0.25)
        {q_H2 = 0;}

    /* Defining H2 uptake Rate */
    /* Uptake Rate: kmol_H2/m3/s */
    r_H2 = c_X * q_H2;

    /* Output source terms: kg_H2/m3/s */
    source = -r_H2 * MW_H2;

    /* Storing check values */
    C_UDMI(cell,thread,24) = q_H2*MW_x*3600;			/* Specific uptake rate: kmol_CO/kg_X/s */
    C_UDMI(cell,thread,25) = source;		            /* Uptake Rate: kg_CO/m3/s */ 

    /* Returning Source Term for */	
    dS[eqn]=0.0;
    return source;
}




DEFINE_SOURCE(CO2_PROD, cell, thread, dS, eqn)
{
    /* Define Variables */
    real source;		                /* Source Term: kg_CO2/m3/s */
    real q_CO2;		                    /* Specific prod rate: kg_CO2/kgx/s */
    real r_CO2;		                    /* Prod Rate: kg_CO/m3/s */
    real c_H2; 		                    /* H2 concentration (liquid): kg/m3 */
    real c_CO; 		                    /* CO concentration (liquid): kg/m3 */
    real q_CO_max, q_H2_max; 		    /* Maximum uptake rate: kg_H2/kg_X/s */
    real K_M_CO, K_M_H2;    		    /* Monod constants: kg_H2/m3 */
    real density_L = C_R(cell,thread);  /* Density: kg/m3 */
    real K_I, K_I_CO, I_CO;             /* Inhibition constants */

    /* Calculating Variables */
    /* Maximum uptake rate: */
    q_CO_max = 1.459/3600;				/* molCO/molx/s */
    q_CO_max = q_CO_max / MW_x;			/* mol CO / gx /h*/
    q_H2_max = 2.565/3600;				/* mol H2/molx/s*/
    q_H2_max = q_H2_max / MW_x;			/* mol H2 / gx /h*/

    /* Monod Constants: mM -> kg_CO/m3 */
    K_M_CO = 0.042*MW_CO/1000;		    /* Monod satutration const: mM -> kg_CO/m3 */
    K_I = 0.246 * MW_CO/1000;		    /* Monod Inhibition const: mM -> kg_CO/m3 */
    K_M_H2 = 0.025*MW_H2/1000;		    /* Monod satutration const: mM -> kg_H2/m3 */

    /* Calculating concentration of CO */
    /* CO concentration (liquid): kg/m3 */
    c_CO = density_L * C_YI(cell,thread,0);
    c_H2 = density_L * C_YI(cell,thread,1);

    /* Calculating inhibition of H2 uptake by CO */
    /* Acetate concentration (liquid): kg/m3 */
    K_I_CO = 0.025 * MW_CO / 1000;			/* Inhibition constant (kg/m3) */
    I_CO = 1 / (1 + c_CO / K_I_CO);			/* Inhibition term (-) */

    /* Calculating uptake Rate */
    /* Specific uptake rate: kmol_H2/kg_X/s */
    q_CO = (q_CO_max * c_CO)/ (c_CO + K_M_CO + (pow(c_CO,2)/K_I));
    q_H2 = (q_H2_max * c_H2)/ (c_H2 + K_M_H2) * I_CO; 

    /* Limiting Biomass Uptake in Droplets and Spray */
    if (C_VOF(cell,thread) <= 0.25)
        {q_H2 = 0;
         q_CO = 0; }


    /* Production rates */
    q_CO2 = 4.0/6.0*q_CO - 1.0/3.0*q_H2;
    r_CO2 = c_X * q_CO2;

    /* Output source terms: kg_CO/m3/s */
    source = r_CO2 * MW_CO2;

    /* Storing check values */
    C_UDMI(cell,thread,26) = q_CO2*MW_x*3600;			/* Specific uptake rate: kmol_CO2/kg_X/s */
    C_UDMI(cell,thread,27) = source;	            	/* Uptake Rate: kg_CO2/m3/s */ 

    /* Returning Source Term for */	
    dS[eqn]=0.0;
return source;
}