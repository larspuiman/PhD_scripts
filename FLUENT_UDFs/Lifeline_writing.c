#include "udf.h"
#include "dpm.h"
#include "dpm_types.h"
#include "mem.h"

/* These scripts is script writes position and velocity of a particle to a file. Executed every timestep in FLUENT. Can easily be augmented with extra variables (pressure, liquid velocity etc) */
DEFINE_EXECUTE_AT_END(particle_data)
{

/* initializations for loop over injection */
Injection *I;
Injection *dpm_injections = Get_dpm_injections();
/* Tracked_Particle *tp;  */
Particle *p; 
FILE *PartOutput;
char PATHNAME_Part[120];				/* provide folder name */
char PATHNAME_Part_full[120];			/* stores path to output directory particle tracks */
int pID;								/* particle ID */
cell_t c;
Thread *t;	
real time, x, y, z;
real clCO, clH2, clCO2;
real qCO, qH2, qCO2;

/* Define global constants */
real MW_CO = 28.01;		/* Molecular mass: g/mol */
real MW_H2 = 2.016;		/* Molecular mass: g/mol */
real MW_CO2 = 44.01;	/* Molecular mass: g/mol */

/* sprintf(PATHNAME_Part,"/tudelft/larspuiman/staff-umbrella/MicroSynC/Lars/PhD/FLUENT/Gaslift_scaledown/DPM/Run2_staticfield/Lifelines");   folder name */
sprintf(PATHNAME_Part,"/home/larspuiman/Documents/Gaslift_scale_down/paper_scale_down/FLUENT/DPM/Run11_25gl_CO2_prod/Lifelines"); /*  folder name */
	loop(I, dpm_injections)																/* Well then, there we go. loop over injections */
		{
			
			loop(p, I->p)
			{

				pID = (int)PP_ID(p);											/* store ID for filename seeking */
				{
					c = PP_CELL(p);
					t = PP_CELL_THREAD(p);
					
					/* load variables */
					time = PP_TIME(p);
				 	x = PP_POS(p)[0];  
					y = PP_POS(p)[1];
					z = PP_POS(p)[2];
					clCO = C_UDMI(c,t,2)/MW_CO*1000;                            /* mol/m3 */
					clH2 = C_UDMI(c,t,9)/MW_H2*1000;                            /* mol/m3 */
					clCO2 = C_UDMI(c,t,14)/MW_CO2*1000;                         /* mol/m3 */
					qCO = C_UDMI(c,t,22);                                       /* mol/molx/h */
					qH2 = C_UDMI(c,t,24);                                       /* mol/molx/h */
					qCO2 = C_UDMI(c,t,26);                                      /* mol/molx/h */

					sprintf(PATHNAME_Part_full,"%s/lifeline_%d.out",PATHNAME_Part,pID);     /* set the filename for this particular particle */
					
    				/* Message("%s \n",PATHNAME_Part_full); */
					/* Message("Particle ID = %d, qCO2 = %f \n", pID, qCO2); */

					PartOutput = fopen(PATHNAME_Part_full,"a");								/* open the file for this particle */
					fprintf(PartOutput,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e \n", time, x, y, z, clCO, clH2, clCO2, qCO, qH2, qCO2);
                    /* fprintf(PartOutput,"%e,%e,%e,%e \n", time, x, y, z); */
					fclose(PartOutput);
    
				
				}

			}

		}
}