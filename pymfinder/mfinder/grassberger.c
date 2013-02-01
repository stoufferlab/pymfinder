/************************************************************************
*
*  File name: grassberger.c
*
*  Description: Functions related to the generation of random networks 
*				using 'go with the winner' algorithm
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
#include "grassberger.h"
#include "random.h"



/*************************** Global variables ****************************/



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
extern int died_colonies;
extern int over_populated_colonies;


/******************************* Functions *******************************/


//clone network entry src to entry dst
int
clone_network(STUBS_RND_NET *SR_ARR,int src,int dst)
{
	int rc=RC_OK;
	int i;

	duplicate_network(SR_ARR[src].N,&SR_ARR[dst].N,"");
	//copy free stubs arrays
	for(i=1;i<=SR_ARR[src].free_sin_stubs;i++){
		SR_ARR[dst].sin_out_stubs[i]=SR_ARR[src].sin_out_stubs[i];
		SR_ARR[dst].sin_in_stubs[i]=SR_ARR[src].sin_in_stubs[i];
	}
	for(i=1;i<=SR_ARR[src].free_dbl_stubs;i++){
		SR_ARR[dst].dbl_out_stubs[i]=SR_ARR[src].dbl_out_stubs[i];
		SR_ARR[dst].dbl_in_stubs[i]=SR_ARR[src].dbl_in_stubs[i];
	}
	SR_ARR[dst].free_sin_stubs=SR_ARR[src].free_sin_stubs;
	SR_ARR[dst].free_dbl_stubs=SR_ARR[src].free_dbl_stubs;
	SR_ARR[dst].active=TRUE;
	return rc;
}



//cloning times are received as input
int 
gen_rand_network_grassberger_method(Network **RN_p, double *clone_time_arr, double *clone_num_p,double *weight_p)
{	
	int rc=RC_OK;
	int s,t;
	int i,j,k;
	STUBS_RND_NET *SR_ARR;
	Network *RN;
	int max_entry=GNRL_ST.r_grass_colony_sz*GNRL_ST.r_grass_colony_max_population;
	int initial_active_entries=GNRL_ST.r_grass_colony_sz*2;
	int out_stub_idx;
	int in_stub_idx;
	int free_stubs;
	int active_net_num;
	int clone_num=0;

	int tot_edges_counter=0;
	int sin_edges_counter=0;
	int dbl_edges_counter=0;

	int self_edge_exists=FALSE;
	//boolien mapping
	int *active_entries;
	int *free_entries;
	int rand_ent;
	int next_clone_time;
	int clone_idx=1;
	double weight,remain_ratio;


	SR_ARR=(STUBS_RND_NET*)calloc(max_entry+1,sizeof(STUBS_RND_NET));
	active_entries=(int*)calloc(max_entry+1,sizeof(int));
	free_entries=(int*)calloc(max_entry+1,sizeof(int));

	for(i=1;i<=max_entry;i++) {
		//allocate
		if(i<=initial_active_entries){
			allocate_network(G_N,&SR_ARR[i].N,"");
		
			SR_ARR[i].active=TRUE;
			active_entries[i]=1;
			free_entries[i]=0;
		}else{
			SR_ARR[i].active=FALSE;
			active_entries[i]=0;
			free_entries[i]=1;
		}

		//gnrl comment:e_dbl_num is twice the number of double edges (they appeasr twice
		SR_ARR[i].free_sin_stubs=G_N->e_sin_num;
		SR_ARR[i].free_dbl_stubs=G_N->e_dbl_num/2;
		
		SR_ARR[i].sin_in_stubs=(int*)calloc(G_N->e_sin_num+1,sizeof(int));
		SR_ARR[i].sin_out_stubs=(int*)calloc(G_N->e_sin_num+1,sizeof(int));
		SR_ARR[i].dbl_in_stubs=(int*)calloc(G_N->e_dbl_num/2+1,sizeof(int));
		SR_ARR[i].dbl_out_stubs=(int*)calloc(G_N->e_dbl_num/2+1,sizeof(int));
		
		//init stubs array
		for(j=1;j<=G_N->e_sin_num;j++){
			SR_ARR[i].sin_out_stubs[j]=G_N->e_arr_sin[j].s;
			SR_ARR[i].sin_in_stubs[j]=G_N->e_arr_sin[j].t;
		}
		for(j=1;j<=G_N->e_dbl_num/2;j++){
			SR_ARR[i].dbl_out_stubs[j]=G_N->e_arr_dbl[2*j].s;
			SR_ARR[i].dbl_in_stubs[j]=G_N->e_arr_dbl[2*j].t;
		}
	}
	next_clone_time=(int)clone_time_arr[1];
	//beginnig the first twice the colony size networks are active
	active_net_num=initial_active_entries;

	
	/* Part 1 - assign double edges */
	
	if(GNRL_ST.out_log_flag && DEBUG_LEVEL==-50){
		fprintf(stdout,"Started assigining Double edges\n");
		fprintf(GNRL_ST.log_fp, "Started assigining Double edges\n");
	}
	free_stubs=SR_ARR[1].free_dbl_stubs;
	while(free_stubs>0){
		tot_edges_counter+=2;
		dbl_edges_counter++;
		for(i=1;i<=max_entry;i++) {
			if(SR_ARR[i].active==FALSE)
				continue;
			//get random out stub and random in stub
			out_stub_idx=get_rand(SR_ARR[i].free_dbl_stubs);
			in_stub_idx=get_rand(SR_ARR[i].free_dbl_stubs);
			s=SR_ARR[i].dbl_out_stubs[out_stub_idx];
			t=SR_ARR[i].dbl_in_stubs[in_stub_idx];
			//check they are not generating conflicts
			//conflicts can be: self edge, existing edge, opposite edge exists
			if((t==s)||(MatGet(SR_ARR[i].N->mat,s,t)>0)||(MatGet(SR_ARR[i].N->mat,t,s)>0)){
				SR_ARR[i].active=FALSE;
				active_entries[i]=0;
				free_entries[i]=1;
				free_network_mem(SR_ARR[i].N);
				free(SR_ARR[i].N);
				active_net_num--;
			}else{
				//add the edge to the building network
				MatAsgn(SR_ARR[i].N->mat,s,t,1);
				MatAsgn(SR_ARR[i].N->mat,t,s,1);

				SR_ARR[i].N->e_arr[tot_edges_counter].s=s;
				SR_ARR[i].N->e_arr[tot_edges_counter].t=t;
				SR_ARR[i].N->e_arr[tot_edges_counter-1].s=t;
				SR_ARR[i].N->e_arr[tot_edges_counter-1].t=s;

				SR_ARR[i].N->e_arr_dbl[(dbl_edges_counter)*2].s=s;
				SR_ARR[i].N->e_arr_dbl[dbl_edges_counter*2].t=t;
				SR_ARR[i].N->e_arr_dbl[dbl_edges_counter*2-1].s=t;
				SR_ARR[i].N->e_arr_dbl[dbl_edges_counter*2-1].t=s;

				//update the stubs arrays
				//move last entry to current
				SR_ARR[i].dbl_out_stubs[out_stub_idx]=SR_ARR[i].dbl_out_stubs[SR_ARR[i].free_dbl_stubs];
				SR_ARR[i].dbl_in_stubs[in_stub_idx]=SR_ARR[i].dbl_in_stubs[SR_ARR[i].free_dbl_stubs];
				SR_ARR[i].free_dbl_stubs--;
			}
		}
		free_stubs--;
		if(GNRL_ST.out_log_flag && DEBUG_LEVEL==-50){
			fprintf(stdout,"%d Edges assigned, Active networks num : %d\n",G_N->e_dbl_num/2-free_stubs,active_net_num);
			fprintf(GNRL_ST.log_fp,"%d Edges assigned, Active networks num : %d\n",G_N->e_dbl_num/2-free_stubs,active_net_num);
		}
		//check how many networks still active
		//if no more active networks free all memory and return RC_COLONY_DIED   
		if(active_net_num==0){
			died_colonies++;
			for(i=1;i<=max_entry;i++) {
				if(SR_ARR[i].active) {
					free_network_mem(SR_ARR[i].N);
					free(SR_ARR[i].N);
				}
				free(SR_ARR[i].sin_in_stubs);
				free(SR_ARR[i].sin_out_stubs);
				free(SR_ARR[i].dbl_in_stubs);
				free(SR_ARR[i].dbl_out_stubs);
			}
			free(SR_ARR);
			free(active_entries);
			free(free_entries);
			printf("\tColony Died\n");
			return (RC_COLONY_DIED);
		}
		//if time to clone
		//here need to check +-1 because any double edge assigned is counted as two edges
		if( (next_clone_time>=tot_edges_counter-1) && (next_clone_time<=tot_edges_counter+1)) {
			if(active_net_num>max_entry/2) {
				over_populated_colonies++;
				for(i=1;i<=max_entry;i++) {
					if(SR_ARR[i].active) {
						free_network_mem(SR_ARR[i].N);
						free(SR_ARR[i].N);
					}
					free(SR_ARR[i].sin_in_stubs);
					free(SR_ARR[i].sin_out_stubs);
					free(SR_ARR[i].dbl_in_stubs);
					free(SR_ARR[i].dbl_out_stubs);
				}
				free(SR_ARR);
				free(active_entries);
				free(free_entries);
				printf("\tOver Populated Colony\n");
				return (RC_COLONY_OVER_POPULATED);
			}
			printf("\tCloning...");
			clone_num++;
			clone_idx++;
			next_clone_time=(int)clone_time_arr[clone_idx];
			if(GNRL_ST.out_log_flag && DEBUG_LEVEL>-50){
				fprintf(stdout,"Cloning networks\n");
				fprintf(GNRL_ST.log_fp,"Cloning networks\n");
				fprintf(stdout,"Active networks before Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(stdout,"%d ",active_entries[j]);
				fprintf(stdout,"\n");
				fprintf(GNRL_ST.log_fp,"Active networks before Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(GNRL_ST.log_fp,"%d ",active_entries[j]);
				fprintf(GNRL_ST.log_fp,"\n");
			}
			for(j=1;j<=max_entry;j++){
				if(active_entries[j]){
					//find next free entry
					k=1;
					while(!free_entries[k]){
						k++;
					}
					//clone network
					clone_network(SR_ARR,j,k);
					free_entries[k]=0;
				}
			}
			//update active entries
			active_net_num=0;
			for(j=1;j<=max_entry;j++){
				if(SR_ARR[j].active==TRUE){
					active_net_num++;
					active_entries[j]=1;
					free_entries[j]=0;
				}
			}
			if(GNRL_ST.out_log_flag && DEBUG_LEVEL>-50){
				fprintf(stdout,"Active networks after Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(stdout,"%d ",active_entries[j]);
				fprintf(stdout,"\n");
				fprintf(GNRL_ST.log_fp,"Active networks after Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(GNRL_ST.log_fp,"%d ",active_entries[j]);
				fprintf(GNRL_ST.log_fp,"\n");
			}
		}
	}

	if(GNRL_ST.out_log_flag && DEBUG_LEVEL==-50){
		fprintf(stdout,"Finished assigining Double edges\n");
		fprintf(GNRL_ST.log_fp, "Finished assigining Double edges\n");
	}

	if(GNRL_ST.out_log_flag && DEBUG_LEVEL==-50){
		fprintf(stdout,"Started assigining Single edges\n");
		fprintf(GNRL_ST.log_fp, "Started assigining Single edges\n");
	}
	

	/* Part 2 - assign single edges */

	free_stubs=SR_ARR[1].free_sin_stubs;
	while(free_stubs>0){
		tot_edges_counter++;
		sin_edges_counter++;
		for(i=1;i<=max_entry;i++) {
			if(SR_ARR[i].active==FALSE)
				continue;
			//get random out stub and random in stub
			out_stub_idx=get_rand(SR_ARR[i].free_sin_stubs);
			in_stub_idx=get_rand(SR_ARR[i].free_sin_stubs);
			s=SR_ARR[i].sin_out_stubs[out_stub_idx];
			t=SR_ARR[i].sin_in_stubs[in_stub_idx];
			//check they are not generating conflicts
			//conflicts can be: self edge, existing edge, opposite edge exists
			if((t==s)||(MatGet(SR_ARR[i].N->mat,s,t)>0)||(MatGet(SR_ARR[i].N->mat,t,s)>0)){
				SR_ARR[i].active=FALSE;
				active_entries[i]=0;
				free_entries[i]=1;
				free_network_mem(SR_ARR[i].N);
				free(SR_ARR[i].N);
				active_net_num--;
			}else{
				//add the edge to the building network
				MatAsgn(SR_ARR[i].N->mat,s,t,1);
				SR_ARR[i].N->e_arr[tot_edges_counter].s=s;
				SR_ARR[i].N->e_arr[tot_edges_counter].t=t;
				
				SR_ARR[i].N->e_arr_sin[sin_edges_counter].s=s;
				SR_ARR[i].N->e_arr_sin[sin_edges_counter].t=t;
				//update the stubs arrays
				//move last entry to current
				SR_ARR[i].sin_out_stubs[out_stub_idx]=SR_ARR[i].sin_out_stubs[SR_ARR[i].free_sin_stubs];
				SR_ARR[i].sin_in_stubs[in_stub_idx]=SR_ARR[i].sin_in_stubs[SR_ARR[i].free_sin_stubs];
				SR_ARR[i].free_sin_stubs--;
			}
		}
		free_stubs--;
		if(GNRL_ST.out_log_flag && DEBUG_LEVEL==-50){
			fprintf(stdout,"%d Edges assigned, Active networks num : %d\n",G_N->e_sin_num-free_stubs,active_net_num);
			fprintf(GNRL_ST.log_fp,"%d Edges assigned, Active networks num : %d\n",G_N->e_sin_num-free_stubs,active_net_num);
		}
		//check how many networks still active
		//if active networks are less then required colony size then failure, return 
		if(active_net_num==0){
		//if(active_net_num<GNRL_ST.r_grass_colony_sz/2){
			//free all mem before returning -1
			died_colonies++;
			for(i=1;i<=max_entry;i++) {
				if(SR_ARR[i].active) {
					free_network_mem(SR_ARR[i].N);
					free(SR_ARR[i].N);
				}
				free(SR_ARR[i].sin_in_stubs);
				free(SR_ARR[i].sin_out_stubs);
				free(SR_ARR[i].dbl_in_stubs);
				free(SR_ARR[i].dbl_out_stubs);
			}
			free(SR_ARR);
			free(active_entries);
			free(free_entries);
			printf("\tColony Died\n");
			return (RC_COLONY_DIED);
		}
		//if time to clone and still there are edges to assign 
		if( (next_clone_time==tot_edges_counter) && (free_stubs>0) ) {
			//if colony over polulated then return right rc 
			if(active_net_num>max_entry/2) {
				over_populated_colonies++;
				for(i=1;i<=max_entry;i++) {
					if(SR_ARR[i].active) {
						free_network_mem(SR_ARR[i].N);
						free(SR_ARR[i].N);
					}
					free(SR_ARR[i].sin_in_stubs);
					free(SR_ARR[i].sin_out_stubs);
					free(SR_ARR[i].dbl_in_stubs);
					free(SR_ARR[i].dbl_out_stubs);
				}
				free(SR_ARR);
				free(active_entries);
				free(free_entries);
				printf("Over Populated Colony\n");
				return (RC_COLONY_OVER_POPULATED);
			}
			printf("\tCloning...");
			clone_num++;
			clone_idx++;
			next_clone_time=(int)clone_time_arr[clone_idx];
			if(GNRL_ST.out_log_flag && DEBUG_LEVEL>-50){
				fprintf(stdout,"Cloning networks\n");
				fprintf(GNRL_ST.log_fp,"Cloning networks\n");
				fprintf(stdout,"Active networks before Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(stdout,"%d ",active_entries[j]);
				fprintf(stdout,"\n");
				for(j=1;j<=max_entry;j++)
					fprintf(stdout,"%d ",SR_ARR[j].active);
				fprintf(stdout,"\n");
				fprintf(GNRL_ST.log_fp,"Active networks before Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(GNRL_ST.log_fp,"%d ",active_entries[j]);
				fprintf(GNRL_ST.log_fp,"\n");
				for(j=1;j<=max_entry;j++)
					fprintf(GNRL_ST.log_fp,"%d ",SR_ARR[j].active);
				fprintf(GNRL_ST.log_fp,"\n");
				
			}
			for(j=1;j<=max_entry;j++){
				if(active_entries[j]){
					k=1;
					while(!free_entries[k]){
						k++;
					}
					//clone network
					clone_network(SR_ARR,j,k);
					free_entries[k]=0;
				}
			}
			//update active entries
			active_net_num=0;
			for(j=1;j<=max_entry;j++){
				if(SR_ARR[j].active==TRUE){
					active_net_num++;
					active_entries[j]=1;
					free_entries[j]=0;
				}
			}
			if(GNRL_ST.out_log_flag && DEBUG_LEVEL>-50){
				fprintf(stdout,"Active networks after Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(stdout,"%d ",active_entries[j]);
				fprintf(stdout,"\n");
				fprintf(GNRL_ST.log_fp,"Active networks after Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(GNRL_ST.log_fp,"%d ",active_entries[j]);
				fprintf(GNRL_ST.log_fp,"\n");
			}
		}
		
	}
	if(GNRL_ST.out_log_flag && DEBUG_LEVEL==-50){
		fprintf(stdout,"Finished assigining Single edges\n");
		fprintf(GNRL_ST.log_fp, "Finished assigining Single edges\n");
	}

	for(i=1;i<=max_entry;i++) {
		if(SR_ARR[i].active) {
			//check there are no Self edges- if there are any then output them to screen and  stop
			for(j=1; j<=G_N->vertices_num; j++) {
				if(MatGet(SR_ARR[i].N->mat,j,j)==1) {
					fprintf(stdout,"Self edges exist in Input Network!!!\n");
					fprintf(stdout,"Self Edges : (%d,%d)\t",j,j);
					self_edge_exists=TRUE;
				}
			}
			if(self_edge_exists==TRUE) {
				printf("Error: Self edges Created in generating random networks\n");
				at_exit(-1);
			}
			
			if((SR_ARR[i].N->edges_num!=(tot_edges_counter)) || (SR_ARR[i].N->e_sin_num!=sin_edges_counter) || (SR_ARR[i].N->e_dbl_num!=2*dbl_edges_counter) ) {
				printf("Degree sequence error in Randomizing network process\n");
				at_exit(-1);
			}
		
		}
	}

	//Now need to choose a random entry from the active entries
	while(!SR_ARR[rand_ent=get_rand(max_entry)].active);
	//duplicate network
	duplicate_network(SR_ARR[rand_ent].N,&RN,"");
	
	//the weight is (2^-c)*remain_ratio
	remain_ratio=(double)active_net_num/(double)initial_active_entries;
	weight=(double)1/((double)(pow(2,clone_num)));
	weight*=remain_ratio;

	//now free all memory
	for(i=1;i<=max_entry;i++) {
		if(SR_ARR[i].active) {
			free_network_mem(SR_ARR[i].N);
			free(SR_ARR[i].N);
		}
		free(SR_ARR[i].sin_in_stubs);
		free(SR_ARR[i].sin_out_stubs);
		free(SR_ARR[i].dbl_in_stubs);
		free(SR_ARR[i].dbl_out_stubs);
		
	}
	free(SR_ARR);
	free(active_entries);
	free(free_entries);

	*RN_p=RN;
	*clone_num_p=clone_num;
	*weight_p=weight;

	return(rc);
}


int 
stats_rand_network_grassberger_method(Network **RN_p, double *clone_time_arr, double *clone_num_p)
{	
	int rc=RC_OK;
	int s,t;
	int i,j,k;
	STUBS_RND_NET *SR_ARR;
	Network *RN;
	int max_entry=GNRL_ST.r_grass_colony_sz*2;
	int out_stub_idx;
	int in_stub_idx;
	int free_stubs;
	int active_net_num;
	int clone_num=0;

	int tot_edges_counter=0;
	int sin_edges_counter=0;
	int dbl_edges_counter=0;

	int self_edge_exists=FALSE;
	//boolien mapping
	int *active_entries;
	int *free_entries;
	int rand_ent;
	int clone_idx=1;

	SR_ARR=(STUBS_RND_NET*)calloc(max_entry+1,sizeof(STUBS_RND_NET));
	active_entries=(int*)calloc(max_entry+1,sizeof(int));
	free_entries=(int*)calloc(max_entry+1,sizeof(int));

	for(i=1;i<=max_entry;i++) {
		//allocate 
		allocate_network(G_N,&SR_ARR[i].N,"");
		SR_ARR[i].active=TRUE;
		active_entries[i]=1;

		//gnrl comment:e_dbl_num is twice the number of double edges (they appeasr twice
		SR_ARR[i].free_sin_stubs=G_N->e_sin_num;
		SR_ARR[i].free_dbl_stubs=G_N->e_dbl_num/2;
		
		SR_ARR[i].sin_in_stubs=(int*)calloc(G_N->e_sin_num+1,sizeof(int));
		SR_ARR[i].sin_out_stubs=(int*)calloc(G_N->e_sin_num+1,sizeof(int));
		SR_ARR[i].dbl_in_stubs=(int*)calloc(G_N->e_dbl_num/2+1,sizeof(int));
		SR_ARR[i].dbl_out_stubs=(int*)calloc(G_N->e_dbl_num/2+1,sizeof(int));
		
		//init stubs array
		for(j=1;j<=G_N->e_sin_num;j++){
			SR_ARR[i].sin_out_stubs[j]=G_N->e_arr_sin[j].s;
			SR_ARR[i].sin_in_stubs[j]=G_N->e_arr_sin[j].t;
		}
		for(j=1;j<=G_N->e_dbl_num/2;j++){
			SR_ARR[i].dbl_out_stubs[j]=G_N->e_arr_dbl[2*j].s;
			SR_ARR[i].dbl_in_stubs[j]=G_N->e_arr_dbl[2*j].t;
		}
	}

	active_net_num=max_entry;

	/* Part 1 - assign double edges */

	if(GNRL_ST.out_log_flag && DEBUG_LEVEL==-50){
		fprintf(stdout,"Started assigining Double edges\n");
		fprintf(GNRL_ST.log_fp, "Started assigining Double edges\n");
	}

	free_stubs=SR_ARR[1].free_dbl_stubs;
	while(free_stubs>0){
		tot_edges_counter+=2;
		dbl_edges_counter++;
		for(i=1;i<=max_entry;i++) {
			if(SR_ARR[i].active==FALSE)
				continue;
			//get random out stub and random in stub
			out_stub_idx=get_rand(SR_ARR[i].free_dbl_stubs);
			in_stub_idx=get_rand(SR_ARR[i].free_dbl_stubs);
			s=SR_ARR[i].dbl_out_stubs[out_stub_idx];
			t=SR_ARR[i].dbl_in_stubs[in_stub_idx];
			//check they are not generating conflicts
			//conflicts can be: self edge, existing edge, opposite edge exists
			if((t==s)||(MatGet(SR_ARR[i].N->mat,s,t)>0)||(MatGet(SR_ARR[i].N->mat,t,s)>0)){
				SR_ARR[i].active=FALSE;
				active_entries[i]=0;
				free_entries[i]=1;
				free_network_mem(SR_ARR[i].N);
				free(SR_ARR[i].N);
				active_net_num--;
			}else{
				//add the edge to the building network
				MatAsgn(SR_ARR[i].N->mat,s,t,1);
				MatAsgn(SR_ARR[i].N->mat,t,s,1);

				SR_ARR[i].N->e_arr[tot_edges_counter].s=s;
				SR_ARR[i].N->e_arr[tot_edges_counter].t=t;
				SR_ARR[i].N->e_arr[tot_edges_counter-1].s=t;
				SR_ARR[i].N->e_arr[tot_edges_counter-1].t=s;

				SR_ARR[i].N->e_arr_dbl[(dbl_edges_counter)*2].s=s;
				SR_ARR[i].N->e_arr_dbl[dbl_edges_counter*2].t=t;
				SR_ARR[i].N->e_arr_dbl[dbl_edges_counter*2-1].s=t;
				SR_ARR[i].N->e_arr_dbl[dbl_edges_counter*2-1].t=s;

				//update the stubs arrays
				//move last entry to current
				SR_ARR[i].dbl_out_stubs[out_stub_idx]=SR_ARR[i].dbl_out_stubs[SR_ARR[i].free_dbl_stubs];
				SR_ARR[i].dbl_in_stubs[in_stub_idx]=SR_ARR[i].dbl_in_stubs[SR_ARR[i].free_dbl_stubs];
				SR_ARR[i].free_dbl_stubs--;
			}
		}
		free_stubs--;
		if(GNRL_ST.out_log_flag && DEBUG_LEVEL==-50){
			fprintf(stdout,"%d Edges assigned, Active networks num : %d\n",G_N->e_dbl_num/2-free_stubs,active_net_num);
			fprintf(GNRL_ST.log_fp,"%d Edges assigned, Active networks num : %d\n",G_N->e_dbl_num/2-free_stubs,active_net_num);
		}
		//check how many networks still active
		//if active networks are less then required colony size then failure, return 
		if(active_net_num==0){
		//if(active_net_num<GNRL_ST.r_grass_colony_sz/2) {
			//free all mem before returning -1
			died_colonies++;
			for(i=1;i<=max_entry;i++) {
				if(SR_ARR[i].active) {
					free_network_mem(SR_ARR[i].N);
					free(SR_ARR[i].N);
				}
				free(SR_ARR[i].sin_in_stubs);
				free(SR_ARR[i].sin_out_stubs);
				free(SR_ARR[i].dbl_in_stubs);
				free(SR_ARR[i].dbl_out_stubs);
			}
			free(SR_ARR);
			free(active_entries);
			free(free_entries);
			return (-1);
		}
		//if less then half of maximum then clone
		if(active_net_num<=max_entry/2) {
			clone_num++;
			//update clone_time_arr
			clone_time_arr[clone_idx++]=tot_edges_counter;
			if(GNRL_ST.out_log_flag && DEBUG_LEVEL>-50){
				fprintf(stdout,"Cloning networks\n");
				fprintf(GNRL_ST.log_fp,"Cloning networks\n");
				fprintf(stdout,"Active networks before Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(stdout,"%d ",active_entries[j]);
				fprintf(stdout,"\n");
				fprintf(GNRL_ST.log_fp,"Active networks before Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(GNRL_ST.log_fp,"%d ",active_entries[j]);
				fprintf(GNRL_ST.log_fp,"\n");
			}
			for(j=1;j<=max_entry;j++){
				if(active_entries[j]){
					//start from middle of the array
					k=1;
					while(!free_entries[k]){
						k++;
					}
					//clone network
					clone_network(SR_ARR,j,k);
					free_entries[k]=0;
				}
			}
			//update active entries
			active_net_num=0;
			for(j=1;j<=max_entry;j++){
				if(SR_ARR[j].active==TRUE){
					active_net_num++;
					active_entries[j]=1;
					free_entries[j]=0;
				}
			}
			if(GNRL_ST.out_log_flag && DEBUG_LEVEL>-50){
				fprintf(stdout,"Active networks after Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(stdout,"%d ",active_entries[j]);
				fprintf(stdout,"\n");
				fprintf(GNRL_ST.log_fp,"Active networks after Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(GNRL_ST.log_fp,"%d ",active_entries[j]);
				fprintf(GNRL_ST.log_fp,"\n");
			}
		}
	}
	if(GNRL_ST.out_log_flag && DEBUG_LEVEL==-50){
		fprintf(stdout,"Finished assigining Double edges\n");
		fprintf(GNRL_ST.log_fp, "Finished assigining Double edges\n");
	}


	/* Part 2 - assign single edges */

	//beginnig all networks  are active
	
	if(GNRL_ST.out_log_flag && DEBUG_LEVEL==-50){
		fprintf(stdout,"Started assigining Single edges\n");
		fprintf(GNRL_ST.log_fp, "Started assigining Single edges\n");
	}
	
	free_stubs=SR_ARR[1].free_sin_stubs;
	while(free_stubs>0){
		tot_edges_counter++;
		sin_edges_counter++;
		for(i=1;i<=max_entry;i++) {
			if(SR_ARR[i].active==FALSE)
				continue;
			//get random out stub and random in stub
			out_stub_idx=get_rand(SR_ARR[i].free_sin_stubs);
			in_stub_idx=get_rand(SR_ARR[i].free_sin_stubs);
			s=SR_ARR[i].sin_out_stubs[out_stub_idx];
			t=SR_ARR[i].sin_in_stubs[in_stub_idx];
			//check they are not generating conflicts
			//conflicts can be: self edge, existing edge, opposite edge exists
			if((t==s)||(MatGet(SR_ARR[i].N->mat,s,t)>0)||(MatGet(SR_ARR[i].N->mat,t,s)>0)){
				SR_ARR[i].active=FALSE;
				active_entries[i]=0;
				free_entries[i]=1;
				free_network_mem(SR_ARR[i].N);
				free(SR_ARR[i].N);
				active_net_num--;
			}else{
				//add the edge to the building network
				MatAsgn(SR_ARR[i].N->mat,s,t,1);
				SR_ARR[i].N->e_arr[tot_edges_counter].s=s;
				SR_ARR[i].N->e_arr[tot_edges_counter].t=t;
				
				SR_ARR[i].N->e_arr_sin[sin_edges_counter].s=s;
				SR_ARR[i].N->e_arr_sin[sin_edges_counter].t=t;
				//update the stubs arrays
				//move last entry to current
				SR_ARR[i].sin_out_stubs[out_stub_idx]=SR_ARR[i].sin_out_stubs[SR_ARR[i].free_sin_stubs];
				SR_ARR[i].sin_in_stubs[in_stub_idx]=SR_ARR[i].sin_in_stubs[SR_ARR[i].free_sin_stubs];
				SR_ARR[i].free_sin_stubs--;
			}
		}
		free_stubs--;
		if(GNRL_ST.out_log_flag && DEBUG_LEVEL==-50){
			fprintf(stdout,"%d Edges assigned, Active networks num : %d\n",G_N->e_sin_num-free_stubs,active_net_num);
			fprintf(GNRL_ST.log_fp,"%d Edges assigned, Active networks num : %d\n",G_N->e_sin_num-free_stubs,active_net_num);
		}
		//check how many networks still active
		//if active networks are less then required colony size then failure, return 
		if(active_net_num==0){
		//if(active_net_num<GNRL_ST.r_grass_colony_sz/2){
			//free all mem before returning -1
			died_colonies++;
			for(i=1;i<=max_entry;i++) {
				if(SR_ARR[i].active) {
					free_network_mem(SR_ARR[i].N);
					free(SR_ARR[i].N);
				}
				free(SR_ARR[i].sin_in_stubs);
				free(SR_ARR[i].sin_out_stubs);
				free(SR_ARR[i].dbl_in_stubs);
				free(SR_ARR[i].dbl_out_stubs);
			}
			free(SR_ARR);
			free(active_entries);
			free(free_entries);
			return (-1);
		}
		//if less then half of maximum then clone
		if(active_net_num<=max_entry/2) {
			clone_num++;
			//update clone_time_arr
			clone_time_arr[clone_idx++]=tot_edges_counter;
			if(GNRL_ST.out_log_flag && DEBUG_LEVEL>-50){
				fprintf(stdout,"Cloning networks\n");
				fprintf(GNRL_ST.log_fp,"Cloning networks\n");
				fprintf(stdout,"Active networks before Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(stdout,"%d ",active_entries[j]);
				fprintf(stdout,"\n");
				for(j=1;j<=max_entry;j++)
					fprintf(stdout,"%d ",SR_ARR[j].active);
				fprintf(stdout,"\n");
				fprintf(GNRL_ST.log_fp,"Active networks before Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(GNRL_ST.log_fp,"%d ",active_entries[j]);
				fprintf(GNRL_ST.log_fp,"\n");
				for(j=1;j<=max_entry;j++)
					fprintf(GNRL_ST.log_fp,"%d ",SR_ARR[j].active);
				fprintf(GNRL_ST.log_fp,"\n");
				
			}
			for(j=1;j<=max_entry;j++){
				if(active_entries[j]){
					//start from middle of the array
					k=1;
					while(!free_entries[k]){
						k++;
					}
					//clone network
					clone_network(SR_ARR,j,k);
					free_entries[k]=0;
				}
			}
			//update active entries
			active_net_num=0;
			for(j=1;j<=max_entry;j++){
				if(SR_ARR[j].active==TRUE){
					active_net_num++;
					active_entries[j]=1;
					free_entries[j]=0;
				}
			}
			if(GNRL_ST.out_log_flag && DEBUG_LEVEL>-50){
				fprintf(stdout,"Active networks after Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(stdout,"%d ",active_entries[j]);
				fprintf(stdout,"\n");
				fprintf(GNRL_ST.log_fp,"Active networks after Cloning\n");
				for(j=1;j<=max_entry;j++)
					fprintf(GNRL_ST.log_fp,"%d ",active_entries[j]);
				fprintf(GNRL_ST.log_fp,"\n");
			}
		}
		
	}
	if(GNRL_ST.out_log_flag && DEBUG_LEVEL==-50){
		fprintf(stdout,"Finished assigining Single edges\n");
		fprintf(GNRL_ST.log_fp, "Finished assigining Single edges\n");
	}


	for(i=1;i<=max_entry;i++) {
		if(SR_ARR[i].active) {
			//check there are no Self edges- if there are any then output them to screen and  stop
			for(j=1; j<=G_N->vertices_num; j++) {
				if(MatGet(SR_ARR[i].N->mat,j,j)==1) {
					fprintf(stdout,"Self edges exist in Input Network!!!\n");
					fprintf(stdout,"Self Edges : (%d,%d)\t",j,j);
					self_edge_exists=TRUE;
				}
			}
			if(self_edge_exists==TRUE) {
				printf("Error: Self edges Created in generating random networks\n");
				at_exit(-1);
			}
			
			if((SR_ARR[i].N->edges_num!=(tot_edges_counter)) || (SR_ARR[i].N->e_sin_num!=sin_edges_counter) || (SR_ARR[i].N->e_dbl_num!=2*dbl_edges_counter) ) {
				printf("Degree sequence error in Randomizing network process\n");
				at_exit(-1);
			}
		
		}
	}

	//Now need to choose a random entry from the active entries
	while(!SR_ARR[rand_ent=get_rand(max_entry)].active);
	//duplicate network
	duplicate_network(SR_ARR[rand_ent].N,&RN,"");
	//now free all memory
	for(i=1;i<=max_entry;i++) {
		if(SR_ARR[i].active) {
			free_network_mem(SR_ARR[i].N);
			free(SR_ARR[i].N);
		}
		free(SR_ARR[i].sin_in_stubs);
		free(SR_ARR[i].sin_out_stubs);
		free(SR_ARR[i].dbl_in_stubs);
		free(SR_ARR[i].dbl_out_stubs);
		
	}
	free(SR_ARR);
	free(active_entries);
	free(free_entries);

	*RN_p=RN;
	*clone_num_p=clone_num;
	return(rc);
}



