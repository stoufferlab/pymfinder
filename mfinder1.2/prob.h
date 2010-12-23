/************************************************************************
*
*  File name: prob.h
*
*  Description: Header file
*  
*  Copyright © 2002-2004 Weizmann Institute of Science, 
*			   76100 Rehovot Israel, All rights reserved 
*
*************************************************************************/
#ifndef __PROB_H
#define __PROB_H



/******************* Definitions ************************/



//convergness mode difference const
#define CONVERGNESS_DIFF_CONST  0.10
//conccentration threshold for taking into acount in the
//convergness mode. subgraphs
//conc in mili (means 0.01 mili)
#define CONC_THRESHOLD          0.01 


/******************* Structures *************************/


/******************* Prototypes *************************/



int
search_subset_prob(Network *N, list *vrtx_set, int mtf_sz, list64 *res, int net_type);
int
count_subgraphs_prob(Network *N, int mtf_sz, list64 **res_p, int net_type, int rand_net_indx);

void
init_edge_pd(Network *N);



#endif
