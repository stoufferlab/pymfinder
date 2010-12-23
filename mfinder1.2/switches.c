/************************************************************************
*
*  File name: switches.c
*
*  Description: Functions related to the generation of random networks 
*				using 'switches' algorithm
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
#include "switches.h"
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



/******************************* Functions *******************************/


/********************************************************
* function : gen_rand_network_switches_method_conserve_double_edges
*	Genreate random network conserving the single node statistics
*   conserving mutual edges
*   Use switch method
*   The original network is taken from G_N structure 
* arguments:
*   RN_p - reference to the random network structure (to be allocated)
*   
* return values:
*	RC_OK - if not error occured
*   RC_ERR - if error occured
*********************************************************/
int
gen_rand_network_switches_method_conserve_double_edges(Network **RN_p, double *switch_ratio)
{
	int num_of_switchs;
	int s1,s2,t1,t2;
	int i,j,k,w;
	int twin_i,twin_j;
	Network *RN;
	int rand_network_checked=FALSE;
	int tot_switches=0;
	int success_switches=0;
	double switches_range;
	
	while (!rand_network_checked) {
		//duplicate source network
		duplicate_network(G_N,&RN,"random_network");
		if (DEBUG_LEVEL>11)
			dump_network(stdout,RN);
		
		//Switch double edges
		
		if(RN->e_dbl_num<=1)
			num_of_switchs=0;
		else {
			// num of switches is round(10*(1+rand)*#edges)
			//num of edges in undirfecetd is actually *2 therefore divided by 2
			switches_range=2*GNRL_ST.r_switch_factor*RN->e_dbl_num/2;
			num_of_switchs=(int)switches_range+get_rand((int)switches_range);
		}
		w=0;//all trials
		k=0;//successfull switches

		if (GNRL_ST.out_log_flag==TRUE)
			fprintf(GNRL_ST.log_fp,"Double edges total number of switches to do : %d\n", num_of_switchs); 
		
		while( (!GNRL_ST.r_global_switch_mode && w<num_of_switchs )){
			w++;
			if( !(w % 50000) && (GNRL_ST.out_log_flag==TRUE) )
				fprintf(GNRL_ST.log_fp,"Double edges Switches: success : %d tries : %d\n",k,w);
			// random edges to choose :
			//only when finding a proper pair of edges make the switch
			i=get_rand(RN->e_dbl_num);
			j=get_rand(RN->e_dbl_num);
			if(i==j)
				continue;
			//this is needed to know if the twin is i-1 or i+1
			if(i & 0x1) 
				twin_i=i+1;
			else
				twin_i=i-1;
			if(j & 0x1)
				twin_j=j+1;
			else
				twin_j=j-1;
			
			s1=RN->e_arr_dbl[i].s;
			t1=RN->e_arr_dbl[i].t;
			s2=RN->e_arr_dbl[j].s;
			t2=RN->e_arr_dbl[j].t;
			
			
			//check : 1.that there are no crossing edges in the network,
			//		  2.there are no common vertices of these 2 edges
			//if hasnt passed the check then have to find other pair to switch
			if ( !( MatGet(RN->mat,s1,t2) || MatGet(RN->mat,s2,t1) || MatGet(RN->mat,t1,s2)|| MatGet(RN->mat,t2,s1)) 
				&& (s1!=t2) &&(s2!=t1) && (s1!=s2) && (t1!=t2) ) {
				if(DEBUG_LEVEL==1 && GNRL_ST.out_log_flag) {
					fprintf(GNRL_ST.log_fp,"switch #%d: (%d %d) (%d %d)\n", k,s1,t1,s2,t2);
					fprintf(GNRL_ST.log_fp,"i: %d, twin i %d , j: %d twin j : %d\n", i,twin_i,j,twin_j);
				}
				//make the switch
				//update edges array
				RN->e_arr_dbl[i].t=t2;
				RN->e_arr_dbl[twin_i].s=t2;
				RN->e_arr_dbl[j].t=t1;
				RN->e_arr_dbl[twin_j].s=t1;
				if(DEBUG_LEVEL==1 && GNRL_ST.out_log_flag) {
					fprintf(GNRL_ST.log_fp,"(%d %d) (%d %d) (%d %d) (%d %d)\n",
						RN->e_arr_dbl[i].s,RN->e_arr_dbl[i].t,
						RN->e_arr_dbl[twin_i].s,RN->e_arr_dbl[twin_i].t,
						RN->e_arr_dbl[j].s,RN->e_arr_dbl[j].t,
						RN->e_arr_dbl[twin_j].s,RN->e_arr_dbl[twin_j].t
						);
				}
				
				//update matrix
				MatAsgn(RN->mat,s1,t1,0);
				MatAsgn(RN->mat,t1,s1,0);
				MatAsgn(RN->mat,s2,t2,0);
				MatAsgn(RN->mat,t2,s2,0);
				if(DEBUG_LEVEL==1 && GNRL_ST.out_log_flag)
					dump_network(GNRL_ST.log_fp,RN);
				MatAsgn(RN->mat,s1,t2,1);
				MatAsgn(RN->mat,t2,s1,1);
				MatAsgn(RN->mat,s2,t1,1);
				MatAsgn(RN->mat,t1,s2,1);
				if(DEBUG_LEVEL==1 && GNRL_ST.out_log_flag)
					dump_network(GNRL_ST.log_fp,RN);
				k++;
			}
		}
		if (GNRL_ST.out_log_flag==TRUE)
			fprintf(GNRL_ST.log_fp,"Double Edges: Finished\n");
		tot_switches+=w;
		success_switches+=k;

		//Switch single edges	
		
		if(RN->e_sin_num<=1)
			num_of_switchs=0;
		else {
			// num of switches is round(10*(1+rand)*#edges)
			switches_range=GNRL_ST.r_switch_factor*RN->e_sin_num;
			num_of_switchs=(int)switches_range+get_rand((int)switches_range);
		}
		w=0; //num of trials 
		k=0; //successfull trials
		if (GNRL_ST.out_log_flag==TRUE)
			fprintf(GNRL_ST.log_fp,"Single edges total number of switches to do : %d\n", num_of_switchs); 
		
		while( (!GNRL_ST.r_global_switch_mode && w<num_of_switchs )){
			w++;
			if(!(w % 50000) && (GNRL_ST.out_log_flag==TRUE))
				fprintf(GNRL_ST.log_fp,"Single edges Switches: success : %d tries : %d\n",k,w);
			//get random indexes for edges
			i=get_rand(RN->e_sin_num);
			j=get_rand(RN->e_sin_num);
			if (i==j)
				continue;
			
			s1=RN->e_arr_sin[i].s;
			t1=RN->e_arr_sin[i].t;
			s2=RN->e_arr_sin[j].s;
			t2=RN->e_arr_sin[j].t;
			
			//check : 1.that there are no crossing edges,
			//		  2.there are no common vertices of these 2 edges
			//		  3.there are no self edges 
			//if hasnt passed the check then have to find other pair to switch
			if ( !( MatGet(RN->mat,s1,t2) || MatGet(RN->mat,s2,t1) || MatGet(RN->mat,t2,s1) || MatGet(RN->mat,t1,s2) )
				&& (s1!=s2) && (s1!=t2) && (t1!=s2) && (t1!=t2) && (s1!=t1) && (s2!=t2)) {
				//make the switch
				//update edges array
				RN->e_arr_sin[i].t=t2;
				RN->e_arr_sin[j].t=t1;
				//update matrix
				MatAsgn(RN->mat,s1,t1,0);
				MatAsgn(RN->mat,s2,t2,0);
				MatAsgn(RN->mat,s1,t2,1);
				MatAsgn(RN->mat,s2,t1,1);
				
				k++;
			}
		}
		if(DEBUG_LEVEL>1)
			dump_network(stdout,RN);
		
		if (GNRL_ST.out_log_flag==TRUE)
				fprintf(GNRL_ST.log_fp,"Single Edges: Finished\n");
		
		tot_switches+=w;
		success_switches+=k;
		
		//check that there are no self edges. If there are any then start again
		rand_network_checked = TRUE;
		if(GNRL_ST.calc_self_edges == FALSE){
			for(i=1;i<=RN->vertices_num;i++) {
				if( MatGet(RN->mat,i,i)==1 ) {
					printf("Self edges still exist building random graph again\n");
					rand_network_checked=FALSE;
					free_network_mem(RN);
					break;
				}
			}
		}
	}

	//merge lists to RN->e_arr
	//
	j=0;
	for(i=1;i<=RN->e_dbl_num;i++) {
		RN->e_arr[++j].s=RN->e_arr_dbl[i].s;
		RN->e_arr[j].t=RN->e_arr_dbl[i].t;
	}
	for(i=1;i<=RN->e_sin_num;i++) {
		RN->e_arr[++j].s=RN->e_arr_sin[i].s;
		RN->e_arr[j].t=RN->e_arr_sin[i].t;
	}
	
	if(w==num_of_switchs*DESPAIR_RATIO && num_of_switchs>0) {
		if(GNRL_ST.out_log_flag)
			fprintf(GNRL_ST.log_fp,"Reached Despair ratio\n");
	}
	if(DEBUG_LEVEL>11)
		dump_network(stdout,RN);
	*RN_p=RN;
	*switch_ratio=(double)success_switches/(double)tot_switches;
	if(success_switches==0)
		printf("Warning: Could not make any Switches in generating random network\n"); 
	return RC_OK;
}







/********************************************************
* function : gen_rand_network_switches_method
*	Genreate random network conserving the single node statistics
*   not conserving mutual edges
*   Use switch method
*   The original network is taken from G_N structure 
* arguments:
*   RN_p - reference to the random network structure (to be allocated)
*   
* return values:
*	RC_OK - if not error occured
*   RC_ERR - if error occured
*********************************************************/
int
gen_rand_network_switches_method(Network **RN_p,double *switch_ratio)
{
	int num_of_switchs;
	int i,j,s1,s2,t1,t2;
	int k,w;
	Network *RN;
	int rand_network_checked=FALSE;
	double switches_range;
	int start_edge,l;
	int tot_switches=0,success_switches=0;
	int edges_set_sz=0;
	
	while (!rand_network_checked) {
		//duplicate source network
		duplicate_network(G_N,&RN,"random network");

		//directed graph
		if ( (GNRL_ST.undirected_flag == FALSE) && (GNRL_ST.r_conserve_layers==FALSE) ) {
			// num of switches is round(10*(1+rand)*#edges)
			//num_of_switchs=(GNRL_ST.r_switch_factor+get_rand(GNRL_ST.r_switch_factor))*RN->edges_num;
			switches_range=2*GNRL_ST.r_switch_factor*RN->edges_num;
			num_of_switchs=(int)switches_range+get_rand((int)switches_range);

			k=0,w=0;
			while( !GNRL_ST.r_global_switch_mode && w<num_of_switchs ) {
				w++;
				//get random indexes for edges
				i=get_rand(RN->edges_num); 
				while ((j=get_rand(RN->edges_num)) == i);
				
				s1=RN->e_arr[i].s;
				t1=RN->e_arr[i].t;
				s2=RN->e_arr[j].s;
				t2=RN->e_arr[j].t;
				
				//check : 1.that there are no crossing edges in the network,
				//		  2.there are no common vertices of these 2 edges
				//if hasnt passed the check then have to find other pair to switch
				if ( !(MatGet(RN->mat,s1,t2) || MatGet(RN->mat,s2,t1))
					&&(s1!=s2) && (s1!=t2) && (t1!=s2) && (t1!=t2) ) {
					//make the switch
					//update edges array
					RN->e_arr[i].t=t2;
					RN->e_arr[j].t=t1;
					//update matrix
					MatAsgn(RN->mat,s1,t1,0);
					MatAsgn(RN->mat,s2,t2,0);
					MatAsgn(RN->mat,s1,t2,1);
					MatAsgn(RN->mat,s2,t1,1);
					
					k++;
				}
				
			}
			//update statistics
			tot_switches+=w;
			success_switches+=k;

		}else if(GNRL_ST.r_conserve_layers==FALSE){
			//undirected graph
			
			// num of switches is round(10*(1+rand)*#edges)
			//num of edges in undirfecetd is actually *2 therefore divided by 2
			num_of_switchs=(10+get_rand(10))*RN->edges_num/2;
			k=0,w=0;
			while(!GNRL_ST.r_global_switch_mode && w<num_of_switchs) {
				w++;
				//get random indexes for edges
				i=get_rand(RN->edges_num/2); 
				while ((j=get_rand(RN->edges_num/2)) == i);
				s1=RN->e_arr[2*i-1].s;
				t1=RN->e_arr[2*i-1].t;
				s2=RN->e_arr[2*j-1].s;
				t2=RN->e_arr[2*j-1].t;
				
				//check : 1.that there are no crossing edges in the network,
				//		  2.there are no common vertices of these 2 edges
				//if hasnt passed the check then have to find other pair to switch
				if ( !( MatGet(RN->mat,s1,t2) || MatGet(RN->mat,s2,t1)) && (s1!=t2) &&(s2!=t1) && (s1!=s2) && (t1!=t2) ) {
					//make the switch
					//update edges array
					RN->e_arr[2*i-1].t=t2;
					RN->e_arr[2*i].s=t2;
					RN->e_arr[2*j-1].t=t1;
					RN->e_arr[2*j].s=t1;
					//update matrix
					MatAsgn(RN->mat,s1,t1,0);
					MatAsgn(RN->mat,t1,s1,0);
					MatAsgn(RN->mat,s2,t2,0);
					MatAsgn(RN->mat,t2,s2,0);
					MatAsgn(RN->mat,s1,t2,1);
					MatAsgn(RN->mat,t2,s1,1);
					MatAsgn(RN->mat,s2,t1,1);
					MatAsgn(RN->mat,t1,s2,1);
					k++;
				}
			}
			//update statistics
			tot_switches+=w;
			success_switches+=k;
		}
		//Conserve layers mode
		else if (GNRL_ST.r_conserve_layers==TRUE){
			

			//start index o fthe edges of this layer in RN->e_arr
			start_edge=0;

			//randomize for each layer separately
			for(l=1;l<=RN->num_of_layers;l++){
				
				//if only a single edge in this layer then break
				if(RN->layer_edges_num[l] <=1)
					continue;
				//decide range of coming edges
				if(RN->layers_range==1){
					edges_set_sz=RN->layer_edges_num[l];
				} else if(RN->layers_range==2){
					//if second layer or last layer then range is one
					if(l==1) 
						edges_set_sz=RN->layer_edges_num[l];
					//edges set is from last two layers
					else
						edges_set_sz=RN->layer_edges_num[l-1]+RN->layer_edges_num[l];
				}

				switches_range=GNRL_ST.r_switch_factor*edges_set_sz;
				num_of_switchs=(int)switches_range+get_rand((int)switches_range);
				
				k=0,w=0;
				while( !GNRL_ST.r_global_switch_mode && w<num_of_switchs ) {
					w++;
					//get random indexes for edges
					i=get_rand(edges_set_sz);
			
					while ( (j=get_rand(edges_set_sz)) == i);
					
					s1=RN->e_arr[i+start_edge].s;
					t1=RN->e_arr[i+start_edge].t;
					s2=RN->e_arr[j+start_edge].s;
					t2=RN->e_arr[j+start_edge].t;
					
					//check : 1.that there are no crossing edges in the network,
					//		  2.there are no common vertices of these 2 edges
					//if hasnt passed the check then have to find other pair to switch
					if ( !(MatGet(RN->mat,s1,t2) || MatGet(RN->mat,s2,t1))
						&&(s1!=s2) && (s1!=t2) && (t1!=s2) && (t1!=t2) ) {
						//make the switch
						//update edges array
						RN->e_arr[i+start_edge].t=t2;
						RN->e_arr[j+start_edge].t=t1;
						//update matrix
						MatAsgn(RN->mat,s1,t1,0);
						MatAsgn(RN->mat,s2,t2,0);
						MatAsgn(RN->mat,s1,t2,1);
						MatAsgn(RN->mat,s2,t1,1);
						
						k++;
					}
				}
				//index of the first node in this layer
				if(RN->layers_range==1)
					start_edge+=RN->layer_edges_num[l];
				else if(RN->layers_range==2){
					if(l>1)
						start_edge+=RN->layer_edges_num[l-1];
				}
				//update statistics
				tot_switches+=w;
				success_switches+=k;
			}
		}
		//check that there are no self edges. If there are any then start again
		rand_network_checked = TRUE;
		if(GNRL_ST.calc_self_edges == FALSE){
			for(i=1;i<=RN->vertices_num;i++) {
				if( MatGet(RN->mat,i,i)==1 ) {
					printf("Self edges still exist building random graph again\n");
					rand_network_checked=FALSE;
					free_network_mem(RN);
					break;
				}
			}
		}
	}
	*switch_ratio=(double)success_switches/(double)tot_switches;
	if(success_switches==0)
		printf("Warning: Could not make any Switches in generating random network\n");
	//dump_network(RN);
	*RN_p=RN;
	return RC_OK;
}



