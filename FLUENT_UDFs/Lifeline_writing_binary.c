#include "udf.h"
#include "dpm.h"
#include "dpm_types.h"
#include "mem.h"

DEFINE_EXECUTE_AT_END(write_part_data_binary)
{
    
    int i;
    /* initializations for loop over injection */
	Injection *I;
	Injection *dpm_injections = Get_dpm_injections();
	Particle *p; 
	int pID;								/* particle ID */

	/* definitions for file writing */
    FILE *PartOutput_conc_bin_tend, *PartOutput_rates;
    char PATHNAME_Part[250];				/* provide folder name */
    char PATHNAME_Part_conc_full[250];			/* stores path to output directory particle tracks */
    char PATHNAME_Part_rates_full[250];			/* stores path to output directory particle tracks */
    sprintf(PATHNAME_Part,"./");   			/* folder name */

    /* Define struct which will store all data in binary */
    /* make all double as that is easy to read later on */
    typedef struct bin_conc_data_struct{
        double pID;                                     /* Particle ID - col 1 */
        double time, x, y, z, clCO, clH2, clCO2;        /* cols 2 - 8 */
        double part_conc_array[num_compounds];          /* cols 9 - 20 */
    }bin_conc_data_struct;

    /* Initialize global structure */
    bin_conc_data_struct bin_conc_dat;
    bin_conc_data_struct* bin_conc_dat_pointer = &bin_conc_dat;

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
			pID = (int)PP_ID(p);											/* store ID for filename seeking */
		
            if (pID < num_particles)
            {

                bin_conc_dat.pID = (double)PP_ID(p);
                
                /* where is our current particle? In a cell, in a thread */
                c = PP_CELL(p);
                t = PP_CELL_THREAD(p);
                
                /* load variables of our current particle into the binary struct */
                bin_conc_dat.time = PP_TIME(p);
                bin_conc_dat.x = PP_POS(p)[0];  
                bin_conc_dat.y = PP_POS(p)[1];
                bin_conc_dat.z = PP_POS(p)[2];
                bin_conc_dat.clCO = C_UDMI(c,t,2);                            /* kg/m3 */
                bin_conc_dat.clH2 = C_UDMI(c,t,9);                            /* kg/m3 */
                bin_conc_dat.clCO2 = C_UDMI(c,t,14);                          /* kg/m3 */

                /* define the current concentrations in the particle */
                for(  i = 0; i < num_compounds; ++i)
                {
                    bin_conc_dat.part_conc_array[i] = P_USER_REAL(p, (i+2));
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

            }               
        }

    }
}