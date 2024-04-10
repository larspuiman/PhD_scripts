/* Function to calculate sum of energy dissipation values, stores in UDMI 19 */

#include "udf.h"
#include "mem.h"
#define ZONE_ID 422             /* check in FLUENT */

DEFINE_ADJUST(diss_sum_calculator, d)
{
	double diss_sum_calculated = 0.0;

        Thread *t = Lookup_Thread(d,ZONE_ID);
        Thread *liq = THREAD_SUB_THREAD(t, 0);
        cell_t c;

	/* Get the value of the thread ID from a user-defined Scheme variable */

	/* thread is only used on compute processes */
        /* loop over cells to calculate the sum of the energy dissipation rate in W [kg m2 s-3]  */
            begin_c_loop_int(c,t)
            {
                diss_sum_calculated += C_D(c,liq)*C_R(c,liq)*C_VOF(c,liq)*C_VOLUME(c,t);    /* m2 / s3 * kg/m3_L * m3_L/m3_R * m3_R */
                /* make sure that right thread pointer is used (liquid or mixture...) */
            }
            end_c_loop_int(c,t)

        /*        Message("diss_sum_calculated %f at node %d \n ",diss_sum_calculated, myid); */

    /* get the total energy dissip rate by summing over the individual nodes */
        diss_sum_calculated = PRF_GRSUM1(diss_sum_calculated);

    /* loop over cells and put number in UDMI */
            begin_c_loop_int (c, t)
            {
                C_UDMI(c, t, 19) = diss_sum_calculated;
            }
            end_c_loop_int (c, t)

}