/************************************************************************
*
*  File name: grassberger.h
*
*  Description: Header file
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#ifndef __GRASS_H
#define __GRASS_H



#include <stdio.h>
#include "random.h"

/******************* Definitions ************************/



#define GEN_NETWORK_GRASSBERGER_DISPAIR 100
//maximum size the colony can be in term of ratio to the original size
#define MAX_COLONY_SZ_RATIO 50

//this defines
#define MAX_CLONES_NUM  1000
//num of colonies for statistics about cloning schedual
#define NUM_OF_COLONIES_4_STATS     10

#define RC_COLONY_DIED              -1
#define RC_COLONY_OVER_POPULATED    -2


/******************* Structures *************************/

#if 0
typedef struct {
	Network *N;
	int active;
	int free_sin_stubs;
	int free_dbl_stubs;
	//arrays of active stubs, the value is the node
	int *sin_in_stubs;//array of active stubs, the value is the node
	int *sin_out_stubs;
	int *dbl_in_stubs;
	int *dbl_out_stubs;
}STUBS_RND_NET;
#endif


/******************* Prototypes *************************/



int
stats_rand_network_grassberger_method(Network **RN_p, double *clone_time_arr, double *clone_num_p);
int
gen_rand_network_grassberger_method(Network **RN_p, double *clone_time_arr, double *clone_num_p, double *weight_p);



#endif
