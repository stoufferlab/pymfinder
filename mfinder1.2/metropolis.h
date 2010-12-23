/************************************************************************
*
*  File name: metropolis.h
*
*  Description: Header file
*  
*  Copyright © 2002-2004 Weizmann Institute of Science, 
*			   76100 Rehovot Israel, All rights reserved 
*
*************************************************************************/
#ifndef __METROP_H
#define __METROP_H

#include <stdio.h>


/******************* Definitions ************************/




#define ARR_NUM 20000
#define TFACTR 0.9   // annealing schedule
#define T0 0.001		// initial temperature
#define ETHRESH 0.005		// energy threshhold
#define ITERATION_F 2	// how many num_of_switches multiples to change


/******************* Structures *************************/

typedef struct{
	int first;
	int second;
	int third;
}triple;


/******************* Prototypes *************************/

void zero_neighbours(int *arr,int *arr_len);
void fill_neighbours(Network *N,int s,int t,int direct,int *arr,int *arr_len);
void fill_13(Network *N,int s1,int t1,int s2,int t2,int sub[14],int add[14]);
void fill_13_dbl(Network *N,int s1,int t1,int s2,int t2,int sub[14],int add[14]);
int is_edge(Network *N,int node1,int node2);
int is_edge2(Network *N,int node1,int node2,int s1,int t1,int s2,int t2);
int is_edge2_dbl(Network *N,int node1,int node2,int s1,int t1,int s2,int t2);
int elm(int state);
void refine(int triples[ARR_NUM][4],int *len,int first,int second,int third);
int equal(int *t1,int *t2,int len);
void met_dump_motifs_res(list64 *res, char *network_name,int vec13[14]);
void met_motifs_search_real(Network *N,Res_tbl *met_res_tbl,int vec[14]);
int metrop(double de, double t);
double energy(int target[14],int current[14]);
void update(int target[14],int minus[14],int plus[14],int fwd);
void output13(int vec[14],FILE *fp);
int gen_rand_network_metrop(Network **RN_p,int real_vec13[14]);


#endif
