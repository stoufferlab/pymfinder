/************************************************************************
*
*  File name: stubs.c
*
*  Description: Functions related to the generation of random networks 
*				using 'stubs' algorithm
*  
*  Copyright © 2002-2004 Weizmann Institute of Science, 
*			   76100 Rehovot Israel, All rights reserved 
*
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "common.h"
#include "list.h"
#include "mat.h"
#include "hash.h"
#include "metropolis.h"
#include "stubs.h"
#include "random.h"
#include "output.h"



/*************************** Global variables ****************************/


int self_edge_violence=0;
int edge_exists_violence=0;
int opposite_edge_exists_violence=0;


/******************************* Externs *********************************/
extern void
init_random_seed();

// return random integer in the interval 1 to max_val
extern int get_rand(int max_val);
extern double get_rand_double();
extern void dump_network(FILE *fp,Network *N);

extern Gnrl_st GNRL_ST;
extern int DEBUG_LEVEL;
extern Network *G_N;



/******************************* Functions *******************************/

int gen_rand_network_stubs_method(Network **RN_p)
{	
	int rc=RC_OK;
	int s,t;
	int i;
	int dispair;
	STUBS_RND_NET *SR;
	int out_stub_idx;
	int in_stub_idx;

	int tot_edges_counter=0;
	int sin_edges_counter=0;
	int dbl_edges_counter=0;

	int self_edge_exists=FALSE;

	//allocate 
	SR=(STUBS_RND_NET*)calloc(1,sizeof(STUBS_RND_NET));
	allocate_network(G_N,&SR->N,"");
	SR->active=TRUE;
	
	//gnrl comment:e_dbl_num is twice the number of double edges (they appeasr twice
	SR->free_sin_stubs=G_N->e_sin_num;
	SR->free_dbl_stubs=G_N->e_dbl_num/2;

	SR->sin_in_stubs=(int*)calloc(G_N->e_sin_num+1,sizeof(int));
	SR->sin_out_stubs=(int*)calloc(G_N->e_sin_num+1,sizeof(int));
	SR->dbl_in_stubs=(int*)calloc(G_N->e_dbl_num/2+1,sizeof(int));
	SR->dbl_out_stubs=(int*)calloc(G_N->e_dbl_num/2+1,sizeof(int));
	
	//init stubs array
	for(i=1;i<=G_N->e_sin_num;i++){
		SR->sin_out_stubs[i]=G_N->e_arr_sin[i].s;
		SR->sin_in_stubs[i]=G_N->e_arr_sin[i].t;
	}
	for(i=1;i<=G_N->e_dbl_num/2;i++){
		SR->dbl_out_stubs[i]=G_N->e_arr_dbl[2*i].s;
		SR->dbl_in_stubs[i]=G_N->e_arr_dbl[2*i].t;
	}

	/* Part 1 - assign single edges */
	dispair=0;
	while(SR->free_sin_stubs>0){
		//get random out stub and random in stub
		out_stub_idx=get_rand(SR->free_sin_stubs);
		in_stub_idx=get_rand(SR->free_sin_stubs);
		s=SR->sin_out_stubs[out_stub_idx];
		t=SR->sin_in_stubs[in_stub_idx];
		//check they are not generating conflicts
		//conflicts can be: self edge, existing edge, opposite edge exists
		if((t==s)||(MatGet(SR->N->mat,s,t)>0)||(MatGet(SR->N->mat,t,s)>0)){
			dispair++;
		}else{
			//add the edge to the building network
			MatAsgn(SR->N->mat,s,t,1);
			SR->N->e_arr[++tot_edges_counter].s=s;
			SR->N->e_arr[tot_edges_counter].t=t;

			SR->N->e_arr_sin[++sin_edges_counter].s=s;
			SR->N->e_arr_sin[sin_edges_counter].t=t;
			//update the stubs arrays
			//move last entry to current
			SR->sin_out_stubs[out_stub_idx]=SR->sin_out_stubs[SR->free_sin_stubs];
			SR->sin_in_stubs[in_stub_idx]=SR->sin_in_stubs[SR->free_sin_stubs];
			SR->free_sin_stubs--;
		}
		if(dispair>ADD_EDGES_DISPAIR_RATIO*G_N->e_sin_num){
			//printf("Reached Edges Despair ratio");
			free_network_mem(SR->N);
			free(SR->sin_in_stubs);
			free(SR->sin_out_stubs);
			free(SR->dbl_in_stubs);
			free(SR->dbl_out_stubs);
			free(SR);
			return (-1);
		}
	}

	/* Part 2 - assign double edges */
	dispair=0;
	while(SR->free_dbl_stubs>0){
		//get random out stub and random in stub
		out_stub_idx=get_rand(SR->free_dbl_stubs);
		in_stub_idx=get_rand(SR->free_dbl_stubs);
		s=SR->dbl_out_stubs[out_stub_idx];
		t=SR->dbl_in_stubs[in_stub_idx];
		//check they are not generating conflicts
		//conflicts can be: self edge, existing edge, opposite edge exists
		if((t==s)||(MatGet(SR->N->mat,s,t)>0)||(MatGet(SR->N->mat,t,s)>0)){
			dispair++;
		}else{
			//add the edge to the building network
			//noe double edge
			MatAsgn(SR->N->mat,s,t,1);
			MatAsgn(SR->N->mat,t,s,1);

			SR->N->e_arr_dbl[(++dbl_edges_counter)*2].s=s;
			SR->N->e_arr_dbl[dbl_edges_counter*2].t=t;
			SR->N->e_arr_dbl[dbl_edges_counter*2-1].s=t;
			SR->N->e_arr_dbl[dbl_edges_counter*2-1].t=s;
			
			SR->N->e_arr[++tot_edges_counter].s=s;
			SR->N->e_arr[tot_edges_counter].t=t;
			SR->N->e_arr[++tot_edges_counter].s=t;
			SR->N->e_arr[tot_edges_counter].t=s;

			//update the stubs arrays
			//move last entry to current
			SR->dbl_out_stubs[out_stub_idx]=SR->dbl_out_stubs[SR->free_dbl_stubs];
			SR->dbl_in_stubs[in_stub_idx]=SR->dbl_in_stubs[SR->free_dbl_stubs];
			SR->free_dbl_stubs--;
		}
		if(dispair>ADD_EDGES_DISPAIR_RATIO*G_N->e_dbl_num){
			//printf("Reached Edges Despair ratio");
			free_network_mem(SR->N);
			free(SR->sin_in_stubs);
			free(SR->sin_out_stubs);
			free(SR->dbl_in_stubs);
			free(SR->dbl_out_stubs);
			free(SR);
			return (-1);
		}
	}

	//check there are no Self edges- if there are any then output them to screen and  stop
	for(i=1; i<=SR->N->vertices_num; i++) {
		if(MatGet(SR->N->mat,i,i)==1) {
			fprintf(stdout,"Self edges exist in Input Network!!!\n");
			fprintf(stdout,"Self Edges : (%d,%d)\t",i,i);
			self_edge_exists=TRUE;
		}
	}
	if(self_edge_exists==TRUE) {
		printf("Error: Self edges Created in generating random networks\n");
		at_exit(-1);
	}
	
	
	if((SR->N->edges_num!=(tot_edges_counter)) || (SR->N->e_sin_num!=sin_edges_counter) || (SR->N->e_dbl_num!=2*dbl_edges_counter) ) {
		printf("Degree sequence error in Randomizing network process\n");
		at_exit(-1);
	}
	
	*RN_p=SR->N;
	free(SR->sin_in_stubs);
	free(SR->sin_out_stubs);
	free(SR->dbl_in_stubs);
	free(SR->dbl_out_stubs);
	free(SR);
	
	return(rc);
}




