/* "simple" script for cell counting, sums over all CPUs */
/* important to define correct cell zone number!! */

#include "udf.h" 
#include "prf.h" 
#define ZONE_ID 58

DEFINE_ON_DEMAND(counting)
{
    #if RP_NODE /* exclude data operations from host process, this stuff should be known to the workers */
    Domain *d = Get_Domain(1);
    Thread *t = Lookup_Thread(d,ZONE_ID);
    cell_t c;
    #endif /* RP_NODE */

    /* should be known to the host + nodes */
    int ncount = 0;
    
    
    #if RP_NODE /* data access only on nodes */
    begin_c_loop_int(c,t)
    {
        ncount += 1;        /* adds one per cell */
    }
    end_c_loop_int(c,t) /* works the same way in serial and parallel */

    /* check the results per node (myid) */
    Message("Number of cells %d at node %d \n ",ncount, myid);
    #endif /* RP_NODE */
   
    /* calculate the total number of cells and bring it to the host */
    ncount = PRF_GISUM1(ncount);
    node_to_host_int_1(ncount);  

    #if RP_HOST
    Message("Total Number of cells %d \n ",ncount);
    #endif 
}