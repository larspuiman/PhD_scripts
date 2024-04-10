/*function to compute gas-liquid mass transfer of carbon monoxide (CO)*/

#include "udf.h"
#include "mem.h"

/* Use correction factor for epsilon? */
int SWITCH_EPS_CORR = 1; /* 0 for not taking it into account */
double f_cor;
double diss_sum;


/* Define global constants */
real MW_CO = 28.01;		/* Molecular mass: g/mol */
real MW_H2 = 2.016;		/* Molecular mass: g/mol */
real MW_CO2 = 44.01;	/* Molecular mass: g/mol */

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
	
/*  Local conditions */
    real eG =  C_VOF(cell,gas);         /* m3G/m3R, local gas hold-up */
	real P = C_P(cell,thread) + 101325; /* Pa, Reactor Pressure*/
	real eps = C_D(cell,liq);		    /* Turbulence energy dissipation, m2/s3 */
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
    double Pmin = 506.7;                  	/* Minimum power input W/m3, Roels Heijnen 1980 */
    double Ptot = Pmin * 565;           	/* Total power input W, liquid-based working volume is 565 m3 */

    if(SWITCH_EPS_CORR == 0)
        {f_cor = 1.0;}
    else
        {
            diss_sum = C_UDMI(cell,thread,19); 
            f_cor = Ptot/diss_sum;          
        }

    k_L_CO_LS =  0.4 * pow(D_L_CO,0.5) * pow((eps*f_cor) / nu_l,0.25); /* m/s */
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
	C_UDMI(cell,thread,0) = P;
	C_UDMI(cell,thread,1) = c_s_CO;
	C_UDMI(cell,thread,2) = c_CO;
	C_UDMI(cell,thread,3) = k_L_CO;
	C_UDMI(cell,thread,4) = m_lg_CO;
	C_UDMI(cell,thread,5) = y_CO;
	C_UDMI(cell,thread,6) = a;
	C_UDMI(cell,thread,7) = k_L_CO * a;
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

/*  Local conditions */
    real eG =  C_VOF(cell,gas);         /* m3G/m3R, local gas hold-up */
	real P = C_P(cell,thread) + 101325; /* Pa, Reactor Pressure*/
	real eps = C_D(cell,liq);		    /* Turbulence energy dissipation, m2/s3 */
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
            diss_sum = C_UDMI(cell,thread,19); 
            f_cor = Ptot/diss_sum;          
        }
	k_L_H2_LS =  0.4 * pow(D_L_H2,0.5) * pow((eps*f_cor) / nu_l,0.25); /* m/s */

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
	C_UDMI(cell,thread,8) = c_s_H2;
	C_UDMI(cell,thread,9) = c_H2;
	C_UDMI(cell,thread,10) = m_lg_H2;
	C_UDMI(cell,thread,11) = y_H2;
	C_UDMI(cell,thread,12) = k_L_H2 * a;
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

/*  Local conditions */
    real eG =  C_VOF(cell,gas);         /* m3G/m3R, local gas hold-up */
	real P = C_P(cell,thread) + 101325; /* Pa, Reactor Pressure*/
	real eps = C_D(cell,liq);		    /* Turbulence energy dissipation, m2/s3 */
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
            diss_sum = C_UDMI(cell,thread,19); 
            f_cor = Ptot/diss_sum;          
        }	
        k_L_CO2_LS =  0.4 * pow(D_L_CO2,0.5) * pow((eps*f_cor) / nu_l,0.25); /* m/s */

/* Calculating MT Coefficient using Higbie relation*/
	real u_slip = C_U(cell,gas) - C_U(cell,liq); /* m/s */
	real v_slip = C_V(cell,gas) - C_V(cell,liq); /* m/s */
	real w_slip = C_W(cell,gas) - C_W(cell,liq); /* m/s */
	real vel_slip = pow(pow(u_slip,2) + pow(v_slip,2) + pow(w_slip,2),0.5); /* m/s */
    k_L_CO2_higbie = pow(4 * D_L_CO2 * vel_slip / d_b / 3.14,0.5); /* m/s */

/* Pick maximum kL value */
    if (k_L_CO2_LS > k_L_CO2_higbie)
        {k_L_CO2 = k_L_CO2_LS;
        C_UDMI(cell,thread,20) = 0; }             /* 0 = Lamonst scott relation for kL */
        else
        {k_L_CO2 = k_L_CO2_higbie;
        C_UDMI(cell,thread,20) = 1; }              /* 1 = Higbie for kL */

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

/* Make correction for gas-hold up volume fraction */
    real xc[ND_ND];
    C_CENTROID(xc, cell, thread);
    real Z = xc[2]; real egcor;
    

    if ((C_VOF(cell, gas) > 0.95) && (Z > 20))
    {
        egcor = 0;
    }
    else
    {
        egcor = C_VOF(cell, gas);
    }

/* Storing check values */
	C_UDMI(cell,thread,13) = c_s_CO2;
	C_UDMI(cell,thread,14) = c_CO2;
	C_UDMI(cell,thread,15) = m_lg_CO2;
	C_UDMI(cell,thread,16) = y_CO2;
	C_UDMI(cell,thread,17) = k_L_CO2 * a;
    C_UDMI(cell,thread,18) = f_cor;
    C_UDMI(cell,thread,21) = egcor;

return (m_lg_CO2);
}