/************************************************************************
*
*  File name: stubs.h
*
*  Description: Header file
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#ifndef __STUBS_H
#define __STUBS_H



#include <stdio.h>
#include "random.h"

/******************* Definitions ************************/



/******************* Structures *************************/


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


/******************* Prototypes *************************/



int gen_rand_network_stubs_method(Network **RN_p);


#endif
