#include "udf.h"
#include "dpm.h"
#include "mem.h"

/* These scripts is script writes position and velocity of a particle to a file. Executed every timestep in FLUENT. Can easily be augmented with extra variables (pressure, liquid velocity etc) */
DEFINE_EXECUTE_AT_END(particle_data)
{

/* initializations for loop over injection */
Injection *I;
Injection *dpm_injections = Get_dpm_injections();
Particle *p;
FILE *PartOutput;
char PATHNAME_Part[120];				/* provide folder name */
char PATHNAME_Part_full[120];			/* stores path to output directory particle tracks */
int pID;								/* particle ID */
cell_t c;
Thread *t;	
real time, x, y, z;
real clCO, clH2, clCO2;
sprintf(PATHNAME_Part,"/tudelft/larspuiman/staff-umbrella/MicroSynC/Lars/PhD/FLUENT/Gaslift_syngas_mix/DPM_3mm_5gl/Part_tracks_2");  /* folder name */
	loop(I, dpm_injections)																/* Well then, there we go. loop over injections */
		{
			
			loop(p, I->p)
			{

				pID = (int)TP_ID(p);											/* store ID for filename seeking */
				{
					c = TP_CELL(p);
					t = TP_CELL_THREAD(p);
					
					/* load variables */
					time = P_TIME(p);
					x = P_POS(p)[0];  
					y = P_POS(p)[1];
					z = P_POS(p)[2];
					clCO = C_UDMI(c,t,2);
					clH2 = C_UDMI(c,t,9);
					clCO2 = C_UDMI(c,t,14);


					sprintf(PATHNAME_Part_full,"%s/lifeline_%d.out",PATHNAME_Part,pID);     /* set the filename for this particular particle */
					
					/* Message("%s \n",PATHNAME_Part_full); */
					
					PartOutput = fopen(PATHNAME_Part_full,"a");								/* open the file for this particle */
					fprintf(PartOutput," %e,%e,%e,%e,%e,%e,%e \n", time, x, y, z, clCO, clH2, clCO2);
					fclose(PartOutput);
				
				
				}

			}

		}
}






