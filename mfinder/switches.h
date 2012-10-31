/************************************************************************
*
*  File name: switches.h
*
*  Description: Header file
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#ifndef __SWITCHES_H
#define __SWITCHES_H



#include <stdio.h>
#include "random.h"

/******************* Definitions ************************/



/******************* Structures *************************/





/******************* Prototypes *************************/

int
gen_rand_network_switches_method_conserve_double_edges(Network **RN_p, double *switch_ratio);
int
gen_rand_network_switches_method(Network **RN_p,double *switch_ratio);

#endif
