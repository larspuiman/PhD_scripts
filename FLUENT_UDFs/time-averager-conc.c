/*function to compute a rolling average in a multi-phase process*/



#include "udf.h"
#include "dpm.h"
#include "mem.h"

/*definitions */
int timer = 1;      /* number of timesteps, initialize at 1 */

/* executes at end of every timestep */
DEFINE_EXECUTE_AT_END(Time_average_std)
{
Domain *d;					/* Load domain */
cell_t cell;				/* will loop over all cells c */
Thread *tr;					/* will loop over all threads tr */
double BU, M2;				/* backup variable for computation */
real cCO, cH2, cCO2;
real qCO, qH2, qCO2;
/* if multiphase specify liquid and gas thread. */ 
Thread *SUBT_Liq;			/* sub-thread for liquid */
Thread *SUBT_Gas;			/* sub-thread for gas */
d = Get_Domain(1);

/* loop over all cells to get average and std data per cell */
  thread_loop_c (tr,d)
    {
		
SUBT_Liq = THREAD_SUB_THREAD(tr,0); /* set liquid thread */
SUBT_Gas = THREAD_SUB_THREAD(tr,1); /* set gas thread */

/* this list can be made arbitrarily long. Make sure to refer to right thread (for epsilon, 
in mixture k-epsilon mode this is the mixture level thread tr, in dispersed mode this is the liquid thread)
make sure to include UDM in simulation! */

       begin_c_loop (cell,tr)
         {
/* velocity averages & stdev, liquid */
BU = C_UDMI(cell,tr,30);
C_UDMI(cell,tr,30) = (C_UDMI(cell,tr,30)*(timer-1) + C_U(cell,SUBT_Liq))/(timer) ;    			/* average liquid_vel u direction */
M2 = C_UDMI(cell,tr,36)+((C_U(cell,SUBT_Liq) - BU)*(C_U(cell,SUBT_Liq)-C_UDMI(cell,tr,30)));		/* standard deviation */
C_UDMI(cell,tr,36) = M2/timer;

BU = C_UDMI(cell,tr,31);
C_UDMI(cell,tr,31) = (C_UDMI(cell,tr,31)*(timer-1) + C_V(cell,SUBT_Liq))/(timer) ;    			/* liquid_vel v */
M2 = C_UDMI(cell,tr,37)+((C_V(cell,SUBT_Liq) - BU)*(C_V(cell,SUBT_Liq)-C_UDMI(cell,tr,31)));
C_UDMI(cell,tr,37) = M2/timer;

BU = C_UDMI(cell,tr,32);
C_UDMI(cell,tr,32) = (C_UDMI(cell,tr,32)*(timer-1) + C_W(cell,SUBT_Liq))/(timer) ;    			/* liquid_vel w */
M2 = C_UDMI(cell,tr,38)+((C_W(cell,SUBT_Liq) - BU)*(C_W(cell,SUBT_Liq)-C_UDMI(cell,tr,32)));
C_UDMI(cell,tr,38) = M2/timer;

/* velocity averages & stdev, gas */
BU = C_UDMI(cell,tr,33);
C_UDMI(cell,tr,33) = (C_UDMI(cell,tr,33)*(timer-1) + C_U(cell,SUBT_Gas))/(timer) ;    			/* average gas_vel u  */
M2 = C_UDMI(cell,tr,39)+((C_U(cell,SUBT_Gas) - BU)*(C_U(cell,SUBT_Gas)-C_UDMI(cell,tr,33)));		/* std dev  */
C_UDMI(cell,tr,39) = M2/timer; 

BU = C_UDMI(cell,tr,34);
C_UDMI(cell,tr,34) = (C_UDMI(cell,tr,34)*(timer-1) + C_V(cell,SUBT_Gas))/(timer) ;    			/* gas vel v */
M2 = C_UDMI(cell,tr,40)+((C_V(cell,SUBT_Gas) - BU)*(C_V(cell,SUBT_Gas)-C_UDMI(cell,tr,34)));
C_UDMI(cell,tr,40) = M2/timer;

BU = C_UDMI(cell,tr,35);
C_UDMI(cell,tr,35) = (C_UDMI(cell,tr,35)*(timer-1) + C_W(cell,SUBT_Gas))/(timer) ;    			/* gas vel w */
M2 = C_UDMI(cell,tr,41)+((C_W(cell,SUBT_Gas) - BU)*(C_W(cell,SUBT_Gas)-C_UDMI(cell,tr,35)));
C_UDMI(cell,tr,41) = M2/timer;

/* turbulence averages & stdev, gas */
BU = C_UDMI(cell,tr,42);
C_UDMI(cell,tr,42) = (C_UDMI(cell,tr,42)*(timer-1) + C_K(cell,SUBT_Liq))/(timer) ;    			/* turb kinetic energy - average */
M2 = C_UDMI(cell,tr,45)+((C_K(cell,SUBT_Liq) - BU)*(C_K(cell,SUBT_Liq)-C_UDMI(cell,tr,42)));	/* turb kinetic energy - stdev */
C_UDMI(cell,tr,45) = M2/timer;

C_UDMI(cell,tr,43) = (C_UDMI(cell,tr,43)*(timer-1) + C_D(cell,SUBT_Liq))/(timer) ;    			/* Turb energy dissipation rate - average */
M2 = C_UDMI(cell,tr,46)+((C_D(cell,SUBT_Liq) - BU)*(C_D(cell,SUBT_Liq)-C_UDMI(cell,tr,43))); 	/* Turb energy dissipation rate - stdev */
C_UDMI(cell,tr,46) = M2/timer;

BU = C_UDMI(cell,tr,44);
C_UDMI(cell,tr,44) = (C_UDMI(cell,tr,44)*(timer-1) + C_VOF(cell,SUBT_Gas))/(timer) ;    		/* Gas fraction */
M2 = C_UDMI(cell,tr,47)+((C_VOF(cell,SUBT_Gas) - BU)*(C_VOF(cell,SUBT_Gas)-C_UDMI(cell,tr,44))); 
C_UDMI(cell,tr,47) = M2/timer;


/* conc_averages */
BU = C_UDMI(cell,tr,48);
cCO = C_UDMI(cell,tr, 2);
C_UDMI(cell,tr,48) = (C_UDMI(cell,tr,48)*(timer-1) + cCO )/(timer) ;    			/* CO conc (g/L) - average */
M2 = C_UDMI(cell,tr,51)+(cCO - BU)*(cCO - C_UDMI(cell,tr,48));	                    /* CO conc (g/L) - stdev */
C_UDMI(cell,tr,51) = M2/timer;

BU = C_UDMI(cell,tr,49);
cH2 = C_UDMI(cell,tr, 9);
C_UDMI(cell,tr,49) = (C_UDMI(cell,tr,49)*(timer-1) + cH2)/(timer) ;    			/* H2 conc (g/L) - average */
M2 = C_UDMI(cell,tr,52)+(cH2 - BU)*(cH2 - C_UDMI(cell,tr,49)); 	            /* H2 conc (g/L) - stdev */
C_UDMI(cell,tr,52) = M2/timer;

BU = C_UDMI(cell,tr,50);
cCO2 = C_UDMI(cell,tr, 14);
C_UDMI(cell,tr,50) = (C_UDMI(cell,tr,50)*(timer-1) + cCO2 )/(timer) ;    		/* CO2 conc (g/L) - average */
M2 = C_UDMI(cell,tr,53)+(( cCO2 - BU)*( cCO2 - C_UDMI(cell,tr,50)));            /* CO2 conc (g/L) - stdev */
C_UDMI(cell,tr,53) = M2/timer;

/* q-averages */
BU = C_UDMI(cell,tr,54);
qCO = C_UDMI(cell,tr, 22);
C_UDMI(cell,tr,54) = (C_UDMI(cell,tr,54)*(timer-1) + qCO )/(timer) ;    		/* CO q rate (mol/mol/h) - average */
M2 = C_UDMI(cell,tr,57)+(qCO - BU)*(qCO - C_UDMI(cell,tr,54));	                /* CO q rate (mol/mol/h) - stdev */
C_UDMI(cell,tr,57) = M2/timer;

BU = C_UDMI(cell,tr,55);
qH2 = C_UDMI(cell,tr, 24);
C_UDMI(cell,tr,55) = (C_UDMI(cell,tr,55)*(timer-1) + qH2)/(timer) ;    			/* H2 q rate (mol/mol/h)  - average */
M2 = C_UDMI(cell,tr,58)+(qH2 - BU)*(qH2 - C_UDMI(cell,tr,55)); 	            /* H2 q rate (mol/mol/h)  - stdev */
C_UDMI(cell,tr,58) = M2/timer;

BU = C_UDMI(cell,tr,56);
qCO2 = C_UDMI(cell,tr, 26);
C_UDMI(cell,tr,56) = (C_UDMI(cell,tr,56)*(timer-1) + qCO2 )/(timer) ;    		/* CO2q rate (mol/mol/h)  - average */
M2 = C_UDMI(cell,tr,59)+(( qCO2 - BU)*( qCO2 - C_UDMI(cell,tr,56)));            /* CO2 q rate (mol/mol/h)  - stdev */
C_UDMI(cell,tr,59) = M2/timer;




         }	
       end_c_loop (cell,tr)
    }
timer = timer+1;	

}


