/************************************************************************
*
*  File name: random.h
*
*  Description: Header file
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#ifndef __RANDOM_H
#define __RANDOM_H


/******************* Definitions ************************/


/******************* Structures *************************/


/*random probability distribution */
typedef struct{
	int size;
	double *cumpd_vec;
}Rnd_pd;


/******************* Prototypes *************************/
void
init_random_seed();
int
get_rand(int max_val);
double
get_rand_double();
void
init_rnd_pd(Rnd_pd **rnd_pd_p,int size,double *pd_vec);
void
free_mem_rnd_pd(Rnd_pd **rnd_pd_p);
int
get_rand_pd(Rnd_pd *rnd_pd);


#endif
