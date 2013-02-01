/************************************************************************
*
*  File name: clustering.h
*
*  Description: Header file
*  
*  Copyright © 2002-2004 Weizmann Institute of Science, 
*			   76100 Rehovot Israel, All rights reserved 
*
*************************************************************************/

#include <stdio.h>


/******************* Definitions ************************/




#define ARR_NUM_C 20000





/******************* Prototypes *************************/

void dump_clustering_series(Network *N,FILE *fp);
int clustering_series(Network *N,double *vec);
void fill_node_neighbours(Network *N,int s,int *arr,int *arr_len);
void dump_clustering_series(Network *N,FILE *fp);
void zero_neighbours_clustering(int *arr,int lim);
double energy_clustering(Network *N,double *rand_cluster_vec);
double sin_metrop_clustering_switches(Network **RN,int nover,int nlimit,double init_t,double *rand_cluster_degree);
double dbl_metrop_clustering_switches(Network **RN,int nover,int nlimit,double init_t,double *rand_cluster_degree);
void update_clustering(Network *N,int s1,int t1,int s2,int t2,double *rand_cluster_degree);
void copy_vec(double *target,double *source,int len);
void update_clustering_single_step(Network *N,int node,int *neighbours,int len,int sign,double *rand_cluster_degree,int node_exclude1,int node_exclude2);
int gen_rand_network_metrop_clustering(Network **RN_p);
