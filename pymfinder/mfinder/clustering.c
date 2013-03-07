/************************************************************************
*
*  File name: clustering.c
*
*  Description: metropolis algorithm implementation
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
#include "clustering.h"
#include "metropolis.h"
#include "results.h"
#include "motif_ids.h"
#include "output.h"
#include "random.h"
#include "switches.h"
#include "stubs.h"

/*************************** Global variables ****************************/


/******************************* Externs *********************************/

extern Gnrl_st GNRL_ST;
extern int DEBUG_LEVEL;
extern Network *G_N;



/******************************* Functions *******************************/
/* fills all neighbours of node i in neighbour_vec and
   returnes their number in K */
/* fills all neighbours of node i in neighbour_vec and
   returnes their number in K */
void fill_node_neighbours(Network *N,int s,int *arr,int *arr_len)
{

	list_item *tmp;
	list_item *tmp2;
	int i;
	int add;

	// add all incoming and outgoing edges from/to node s
	tmp = list_get_next(N->mat->spr->m[s].to,NULL);
	tmp2 = list_get_next(N->mat->spr->m[s].from,NULL);

	while (tmp != NULL)
	{
		arr[*arr_len]=tmp->val;
		(*arr_len)=(*arr_len)+1;
		tmp=list_get_next(N->mat->spr->m[s].to,tmp);
	}

	add=1;
	while (tmp2 != NULL)
	{
		for (i=0;i<=*arr_len;i++)
			if (arr[i]==tmp2->val)
				add=0;
		if (add)
		{
			arr[*arr_len]=tmp2->val;
			(*arr_len)=(*arr_len)+1;
		}
		tmp2=list_get_next(N->mat->spr->m[s].from,tmp2);
	}

}

/* output clustering series to a file */
/* output clustering series to a file */


int clustering_series(Network *N,double *vec)
{
	register int i,j,m;
	int neighbour_vec[ARR_NUM_C];
	int count,K,max_K;
	double pairs;
	double C;

	max_K=0;
	for (i=1; i<=N->vertices_num;i++)
	{
		K=0;
		zero_neighbours_clustering(neighbour_vec,max_K);
		fill_node_neighbours(N,i,neighbour_vec,&K);
		if (K>max_K)
			max_K=K;
		count=0;
		// check all pairs of nodes, increment count if connected
		for (j=0;j<max_K;j++)
		{
			for (m=0;m<K;m++)
			{
				if ((is_edge(N,neighbour_vec[j],neighbour_vec[m]))||(is_edge(N,neighbour_vec[m],neighbour_vec[j])))
					count++;
			}
		}
		// finally convert to clustering coefficient
		if (K<2)
			C=-1;
		else
		{
			pairs=(double)(K*(K-1)/2);
			C=((double)(count/2)/(double)(pairs));
		}
		if (vec==NULL)
			N->cluster_degree[i]=C;
		else
			vec[i]=C;

		//printf("node=%d C=%lf\n",i,C);
	}


	return 0;
}




/* output clustering series to a file */

void dump_clustering_series(Network *N,FILE *fp)
{
	register int i;

	fprintf(fp,"\n");
	for (i=1; i<=N->vertices_num;i++)
		fprintf(fp,"%lf\n",N->cluster_degree[i]);

}

void zero_neighbours_clustering(int *arr,int lim)
{
	register int i;

	for (i=0;i<lim+5;i++)
	{
		arr[i]=0;
	}

}

/*  generate random network which preserves clustering series of real network       */

int gen_rand_network_metrop_clustering(Network **RN_p)
{
	int rc=RC_OK;
	int sin_num_of_switchs,dbl_num_of_switchs,num_at_t;
	//int s2,t1,t2;
	int i,j,l,ii;
	int k;  // number of succesful switches
	int w;  // total number of switches
	Network *RN;
	int rand_network_checked=FALSE;
	double E1=1,E2=1;
	int	snover,dnover; // maximum # of graph changes at any temperature
	//int nlimit,snlimit,dnlimit; // maximum # of succesful graph changes at any temperature
	int snlimit,dnlimit; // maximum # of succesful graph changes at any temperature

	double t;  // the temperature
	int success;
	int dispair=0;
	double dummy;
	int switches_range;
	double *switch_ratio_arr;
	double *rand_cluster_vec;

	switch_ratio_arr=(double*)calloc(GNRL_ST.rnd_net_num+1,sizeof(double));

	// rand_cluster_vec will hold clustering series of the random network
	rand_cluster_vec=(double*)calloc(G_N->vertices_num+5,sizeof(double));
	for (i=1;i<=G_N->vertices_num+2;i++)
	{
		rand_cluster_vec[i]=G_N->cluster_degree[i];
		//printf("C=%lf\n",rand_cluster_vec[i]);
	}

	sin_num_of_switchs=0;
	if(G_N->e_sin_num<=1)
			sin_num_of_switchs=0;
		else {
			// num of switches is round(10*(1+rand)*#edges)
			//num of edges in undirfecetd is actually *2 therefore divided by 2
			switches_range=(int)(2*GNRL_ST.r_switch_factor*G_N->e_sin_num);
			sin_num_of_switchs=(int)switches_range+get_rand((int)switches_range);
		}
	if(G_N->e_dbl_num<=1)
			dbl_num_of_switchs=0;
		else {
			// num of switches is round(10*(1+rand)*#edges)
			//num of edges in undirfecetd is actually *2 therefore divided by 2
			switches_range=(int)(2*GNRL_ST.r_switch_factor*G_N->e_dbl_num/2);
			dbl_num_of_switchs=(int)switches_range+get_rand((int)switches_range);
		}

	t=GNRL_ST.t_init;

	if (GNRL_ST.use_stubs_method==TRUE)
	{

		success=0;
		dispair=0;
		while ((success==FALSE)&&(dispair<GEN_NETWORK_DISPAIR))
		{
			//printf("trying ron\n");
			rc=gen_rand_network_stubs_method(&RN);
			if(rc==RC_OK)
				success=TRUE;
			else
				success=FALSE;
			if (success==FALSE)
			{
				//printf("dispair - discard\n");
				dispair++;
			}

		}

		if (dispair==GEN_NETWORK_DISPAIR) // give up try switch method
			gen_rand_network_switches_method_conserve_double_edges(&RN,&dummy);
		else
		{
			success=update_network(RN,G_N);
			if (!success)
			{
				printf("Error - degrees don't match\n");
				at_exit(-1);
			}

			//shuffle_double_edges(RN);

			if (DEBUG_LEVEL>20 && GNRL_ST.out_log_flag)
			{
				fprintf(GNRL_ST.log_fp,"The network :\n\n");
				for (l=1;l<=(*RN).edges_num;l++)
					fprintf(GNRL_ST.log_fp,"%d %d %d\n",(*RN).e_arr[l].t,(*RN).e_arr[l].s,1);
			}
		}
	}

	else
	{
		gen_rand_network_switches_method_conserve_double_edges(&RN, &dummy);
		//dump_network(stdout,RN);
		clustering_series(RN,rand_cluster_vec);
		if(GNRL_ST.out_log_flag)
		{
			dump_network(GNRL_ST.log_fp,RN);
			fflush(GNRL_ST.log_fp);
			for (ii=1; ii<=RN->vertices_num;ii++)
				fprintf(GNRL_ST.log_fp,"C(%d)=%lf\n",ii,RN->cluster_degree[i]);
			fflush(GNRL_ST.log_fp);
		}

	}

	/**************************************************************/
	/**************************************************************/
	/**************************************************************/
	/* Part II - metropolis clustering algorithm. Perform switches with probability*/
	/**************************************************************/
	/**************************************************************/
	/**************************************************************/

	/**************************************************************/
		/* metropolis parameters depend on num_of _switches*/


		E1=energy_clustering(G_N,rand_cluster_vec);

	w=0;

	if (dbl_num_of_switchs==0)
	{
		snover=sin_num_of_switchs*(int)GNRL_ST.iteration_factor;
		snlimit=snover/50;
		t=sin_metrop_clustering_switches(&RN,snover,snlimit,GNRL_ST.t_init,rand_cluster_vec);

	}
	else
	{
		dnover=dbl_num_of_switchs*(int)GNRL_ST.iteration_factor/30;
		dnlimit=dnover/50;
		if (dnover<10)
		{

			dnover=dbl_num_of_switchs;
			dnlimit=dnover;
		}
		t=dbl_metrop_clustering_switches(&RN,dnover,dnlimit,GNRL_ST.t_init,rand_cluster_vec);

		snover=sin_num_of_switchs*(int)GNRL_ST.iteration_factor/30;
		snlimit=snover/50;
		k=0;
		num_at_t=0;
		t=sin_metrop_clustering_switches(&RN,snover,snlimit,GNRL_ST.t_init,rand_cluster_vec);


		if(GNRL_ST.out_metrop_mat_flag) {
			fprintf(GNRL_ST.mat_metrop_fp,"End E\n");
			printf("End E\n");
		}
		for (i=1;i<20;i++)
		{
			E2=energy_clustering(G_N,rand_cluster_vec);
			if (E2>0)
			{
				t=dbl_metrop_clustering_switches(&RN,dnover,dnlimit,GNRL_ST.t_init,rand_cluster_vec);
				t=sin_metrop_clustering_switches(&RN,snover,snlimit,GNRL_ST.t_init,rand_cluster_vec);
			}
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

	if(DEBUG_LEVEL>1 && GNRL_ST.out_metrop_mat_flag) {
		fprintf(GNRL_ST.mat_metrop_fp,"The network :\n");
		for (i=1;i<=RN->edges_num;i++)
			fprintf(GNRL_ST.mat_metrop_fp,"%d %d %d\n",RN->e_arr[i].t,RN->e_arr[i].s,1);
		fprintf(GNRL_ST.mat_metrop_fp,"End of network\n");
		fprintf(GNRL_ST.mat_metrop_fp,"\n\n");
	}
	if(GNRL_ST.out_metrop_mat_flag) {
		fprintf(GNRL_ST.mat_metrop_fp,"Real triadic census :\n ");
//		output13(real_vec13,GNRL_ST.mat_metrop_fp);
		fprintf(GNRL_ST.mat_metrop_fp,"Final random triadic census :\n ");
//		output13(rand_vec13,GNRL_ST.mat_metrop_fp);
	}

	clustering_series(RN,rand_cluster_vec);

	if(GNRL_ST.out_clustering==TRUE)
		dump_clustering_series(G_N,GNRL_ST.clust_fp);

	*RN_p=RN;

	if(GNRL_ST.out_log_flag)
	{
		dump_network(GNRL_ST.log_fp,RN);
		fflush(GNRL_ST.log_fp);
		for (ii=1; ii<=RN->vertices_num;ii++)
			fprintf(GNRL_ST.log_fp,"C(%d)=%lf\n",ii,RN->cluster_degree[ii]);
		fprintf(GNRL_ST.log_fp,"change, E1=%lf, E2=%lf\n",E1,E2);
		fflush(GNRL_ST.log_fp);
	}


	if(rand_cluster_vec!=NULL)
			free(rand_cluster_vec);
	if(switch_ratio_arr!=NULL)
			free(switch_ratio_arr);

	return RC_OK;
}


/* metropolis clustering energy - quantifies distance between N->cluster_degree and the
   random network cluster degree */

double energy_clustering(Network *N,double *rand_cluster_vec)
{

	register int i;
	double result=0.0,t;
	//double mone,mechane,sm,df;
	double sm,df;

	for (i=1;i<=N->vertices_num;i++)
	{
		sm=(double)(N->cluster_degree[i]+rand_cluster_vec[i]);
		if (sm!=0)
		{
			df=(double)(fabs(N->cluster_degree[i]-rand_cluster_vec[i]));
			t=(double)((double)df/(double)sm);
			result+=t;
		}
	}
	return(result);

}


/* perform metropolis iteration exchanging single edges. Output the temperature */
double sin_metrop_clustering_switches(Network **RN,int nover,int nlimit,double init_t,double *rand_cluster_degree)
{

	double t=init_t;
	register int w,k,num_at_t;
	int s1,t1,s2,t2;
	int i,j,l;
	double E1=0;
	double E2=0;
	double dE=0;
	int make_change;
	double *last_cluster_degree;
	if(GNRL_ST.out_metrop_mat_flag==TRUE)
		printf("\nBeginning single switches\n");

	k=0,w=0,num_at_t=0;
	t=init_t;


	last_cluster_degree=(double*)calloc(G_N->vertices_num+5,sizeof(double));

	//printf("starting single swithces\n");
	E1=energy_clustering(*RN,rand_cluster_degree);

	if ((w%1000)==0 && GNRL_ST.out_metrop_mat_flag)
		fprintf(GNRL_ST.mat_metrop_fp,"\nEsin :\n");

	//printf("T=%lf, E=%lf\n",t,E1);
	while ((w<nover)&&(E1>GNRL_ST.e_thresh))
	//while (E1>GNRL_ST.e_thresh)
	{
		w++;
		num_at_t++;

		if ((k>nlimit)||(num_at_t>nlimit))
		{
			num_at_t=0;
			k=0;
			t *= TFACTR;		// lower temperature
			//printf("T=%lf, E=%lf\n",t,E1);
		}

		//get random indexes for edges
		i=get_rand((*RN)->e_sin_num);
		while ((j=get_rand((*RN)->e_sin_num)) == i);

		s1=(*RN)->e_arr_sin[i].s;
		t1=(*RN)->e_arr_sin[i].t;
		s2=(*RN)->e_arr_sin[j].s;
		t2=(*RN)->e_arr_sin[j].t;


		//check : 1.that there are no crossing edges,
		//		  2.there are no common vertices of these 2 edges
		//if hasnt passed the check then have to find other pair to switch

		if ( !(MatGet((*RN)->mat,s1,t2) || MatGet((*RN)->mat,s2,t1)|| MatGet((*RN)->mat,t2,s1) || MatGet((*RN)->mat,t1,s2))
				&&(s1!=s2) && (s1!=t2) && (t1!=s2) && (t1!=t2) )
		{
			copy_vec(last_cluster_degree,rand_cluster_degree,G_N->vertices_num);
			/* make the relevant changes to rand_cluster_degree*/
			update_clustering((*RN),s1,t1,s2,t2,rand_cluster_degree);
			E2=energy_clustering(*RN,rand_cluster_degree);
			dE=E2-E1;
			//	printf("k=%d E=%lf\n",k,E1);
			make_change=metrop(dE,t);


			if (make_change)
			{
				if ((w%1000)==0 && GNRL_ST.out_metrop_mat_flag)
				{
					fprintf(GNRL_ST.mat_metrop_fp,"%lf\n",E1);
					printf("success=%d k=%d E1=%lf T=%lf\n",k,w,E1,t);
				}

				if(DEBUG_LEVEL>1 && GNRL_ST.out_metrop_mat_flag)
				{
					fprintf(GNRL_ST.mat_metrop_fp,"The network :\n\n");
					for (l=1;l<=(*RN)->e_sin_num;l++)
						fprintf(GNRL_ST.mat_metrop_fp,"%d %d %d\n",(*RN)->e_arr_sin[l].t,(*RN)->e_arr_sin[l].s,1);
				}


				E1=E2;
				//make the switch
				//update edges array
				(*RN)->e_arr_sin[i].t=t2;
				(*RN)->e_arr_sin[j].t=t1;
				//update matrix
				MatAsgn((*RN)->mat,s1,t1,0);
				MatAsgn((*RN)->mat,s2,t2,0);
				MatAsgn((*RN)->mat,s1,t2,1);
				MatAsgn((*RN)->mat,s2,t1,1);
				//fprintf(GNRL_ST.mat_metrop_fp,"%lf %lf\n  ",E1,t);
				k++;

				//printf("change, E1=%lf, E2=%lf\n",E1,E2);
			}
			else
			/* change cluster degree vector back */
			{
				copy_vec(rand_cluster_degree,last_cluster_degree,G_N->vertices_num);
				//printf("no change, E1=%lf, E2=%lf\n",E1,E2);
			}
		}
	}
	if(GNRL_ST.out_metrop_mat_flag==TRUE) {
		fprintf(GNRL_ST.mat_metrop_fp,"%lf\n",E1);
		printf("success=%d k=%d E1=%lf T=%lf\n",k,w,E1,t);
	}
	if(last_cluster_degree!=NULL)
			free(last_cluster_degree);
	return(t);
}


/* perform metropolis iteration exchanging double edges. Output the temperature */
double dbl_metrop_clustering_switches(Network **RN,int nover,int nlimit,double init_t,double *rand_cluster_degree)
{

	double t=init_t;
	register int w,k,num_at_t;
	int s1,t1,s2,t2;
	int i,j;
	double E1,E2,dE;
	int make_change;
	int tries; //number of tries to find proper pair of edges to exchange
	int twin_i=0;
	int twin_j=0;
	double *last_cluster_degree;

	last_cluster_degree=(double*)calloc(G_N->vertices_num+5,sizeof(double));

	w=0;
	if ((w%1000)==0 && GNRL_ST.out_metrop_mat_flag)
		fprintf(GNRL_ST.mat_metrop_fp,"\nBeginning double switches\n");
	k=0,w=0,num_at_t=0;
	t=init_t;

	printf("starting mutual swithces\n");
	E1=energy_clustering(*RN,rand_cluster_degree);

	while ((w<nover)&&(E1>GNRL_ST.e_thresh))
	{
		w++;
		num_at_t++;
		if ((k>nlimit)||(num_at_t>nlimit))
		{
			num_at_t=0;
			k=0;
			t *= TFACTR;		// lower temperature
			//printf("T=%lf, E=%lf\n",t,E1);
		}
		w++;
		// random edges to choose : until finding a proper pair of edges which can be candidates
		// should not choose the same edge or the twin one (of the double)
		// try this for 100 times if not getting good candidate then probably
		// there is a only one double edge in the netwprk
		for(i=get_rand((*RN)->e_dbl_num),j=get_rand((*RN)->e_dbl_num),tries=0;
		((j==i) || (j==twin_i)) && (tries<100) ;i=get_rand((*RN)->e_dbl_num),j=get_rand((*RN)->e_dbl_num),tries++);
		//the way out of deadlock
		if(tries==100 || (*RN)->e_dbl_num==2)
			break;
		//this is needed to know if the twin is i-1 or i+1
		if(i & 0x1)
			twin_i=i+1;
		else
			twin_i=i-1;
		if(j & 0x1)
			twin_j=j+1;
		else
			twin_j=j-1;

		s1=(*RN)->e_arr_dbl[i].s;
		t1=(*RN)->e_arr_dbl[i].t;
		s2=(*RN)->e_arr_dbl[j].s;
		t2=(*RN)->e_arr_dbl[j].t;

		//check : 1.that there are no crossing edges in the network,
		//		  2.there are no common vertices of these 2 edges
		//if hasnt passed the check then have to find other pair to switch
		if (( !( MatGet((*RN)->mat,s1,t2) || MatGet((*RN)->mat,s2,t1) || MatGet((*RN)->mat,t1,s2)|| MatGet((*RN)->mat,t2,s1))
			&& (s1!=t2) &&(s2!=t1) && (s1!=s2) && (t1!=t2)) )
		{
			copy_vec(last_cluster_degree,rand_cluster_degree,G_N->vertices_num);
			/* make the relevant changes to rand_cluster_degree*/
			update_clustering((*RN),s1,t1,s2,t2,rand_cluster_degree);
			E2=energy_clustering(*RN,rand_cluster_degree);
			dE=E2-E1;
			//	printf("k=%d E=%lf\n",k,E1);
			make_change=metrop(dE,t);

			if (make_change)
			{
				if(DEBUG_LEVEL==1 && GNRL_ST.out_log_flag) {
					fprintf(GNRL_ST.log_fp,"switch #%d: (%d %d) (%d %d)\n", k,s1,t1,s2,t2);
					fprintf(GNRL_ST.log_fp,"i: %d, twin i %d , j: %d twin j : %d\n", i,twin_i,j,twin_j);
				}
				if ((w%1000)==0 && GNRL_ST.out_metrop_mat_flag)
				{
					fprintf(GNRL_ST.mat_metrop_fp,"%lf\n",E1);
					printf("success=%d k=%d E1=%lf T=%lf\n",k,w,E1,t);
				}
				E1=E2;
				if ((w%1000)==0 && GNRL_ST.out_metrop_mat_flag)
						printf("success=%d k=%d E1=%lf T=%lf\n",k,w,E1,t);
				//make the switch
				//update edges array
				(*RN)->e_arr_dbl[i].t=t2;
				(*RN)->e_arr_dbl[twin_i].s=t2;
				(*RN)->e_arr_dbl[j].t=t1;
				(*RN)->e_arr_dbl[twin_j].s=t1;
				if(DEBUG_LEVEL==1 && GNRL_ST.out_log_flag) {
					fprintf(GNRL_ST.log_fp,"(%d %d) (%d %d) (%d %d) (%d %d)\n",
					(*RN)->e_arr_dbl[i].s,(*RN)->e_arr_dbl[i].t,
					(*RN)->e_arr_dbl[twin_i].s,(*RN)->e_arr_dbl[twin_i].t,
					(*RN)->e_arr_dbl[j].s,(*RN)->e_arr_dbl[j].t,
					(*RN)->e_arr_dbl[twin_j].s,(*RN)->e_arr_dbl[twin_j].t
					);
				}

				//update matrix
				MatAsgn((*RN)->mat,s1,t1,0);
				MatAsgn((*RN)->mat,t1,s1,0);
				MatAsgn((*RN)->mat,s2,t2,0);
				MatAsgn((*RN)->mat,t2,s2,0);
				if(DEBUG_LEVEL==1 && GNRL_ST.out_log_flag)
					dump_network(GNRL_ST.log_fp,(*RN));
				MatAsgn((*RN)->mat,s1,t2,1);
				MatAsgn((*RN)->mat,t2,s1,1);
				MatAsgn((*RN)->mat,s2,t1,1);
				MatAsgn((*RN)->mat,t1,s2,1);
				if(DEBUG_LEVEL==1 && GNRL_ST.out_log_flag)
					dump_network(GNRL_ST.log_fp,(*RN));
				if (DEBUG_LEVEL>1 && GNRL_ST.out_metrop_mat_flag)
				{
					fprintf(GNRL_ST.mat_metrop_fp,"%lf %lf  ",E1,t);
//					output13(rand_vec13,GNRL_ST.mat_metrop_fp);
				}
				k++;
				if(GNRL_ST.out_log_flag)
					dump_network(GNRL_ST.log_fp,(*RN));

			}
			else
			/* change cluster degree vector back */
			{
				copy_vec(rand_cluster_degree,last_cluster_degree,G_N->vertices_num);
			}
		}
	}
	if(GNRL_ST.out_metrop_mat_flag==TRUE){
		fprintf(GNRL_ST.mat_metrop_fp,"%lf\n",E1);
		printf("success=%d k=%d E1=%lf T=%lf\n",k,w,E1,t);
	}
	if(last_cluster_degree!=NULL)
			free(last_cluster_degree);
	return(t);
}




/****************************************************************/
/* update the relevant indices in rand_cluster_degree. These are:
			   1. nodes connected to s1 and t1
			   2. nodes connected to s1 and t2
			   3. nodes connected to s2 and t2
			   4. nodes connected to s2 and t1
			   5. nodes connected to t1 and s1
			   6. nodes connected to t1 and s2
			   7. nodes connected to t2 and s2
			   8. nodes connected to t2 and s1

***************************************************************/
void update_clustering(Network *N,int s1,int t1,int s2,int t2,double *rand_cluster_degree)
{
	int *neighbour_vec_s1;
	int *neighbour_vec_t1;
	int *neighbour_vec_s2;
	int *neighbour_vec_t2;
	int len_n_s1=0;
	int len_n_t1=0;
	int len_n_s2=0;
	int len_n_t2=0;

	neighbour_vec_s1=(int*)calloc(N->vertices_num+5,sizeof(int));
	neighbour_vec_t1=(int*)calloc(N->vertices_num+5,sizeof(int));
	neighbour_vec_s2=(int*)calloc(N->vertices_num+5,sizeof(int));
	neighbour_vec_t2=(int*)calloc(N->vertices_num+5,sizeof(int));


	fill_node_neighbours(N,s1,neighbour_vec_s1,&len_n_s1);
	fill_node_neighbours(N,t1,neighbour_vec_t1,&len_n_t1);
	fill_node_neighbours(N,s2,neighbour_vec_s2,&len_n_s2);
	fill_node_neighbours(N,t2,neighbour_vec_t2,&len_n_t2);

	/* decrement by 1/(k 2) rand_cluster_degree for nodes of classes 1,3,5,7*/
	update_clustering_single_step(N,t1,neighbour_vec_s1,len_n_s1,-1,rand_cluster_degree,-1,-1);
	update_clustering_single_step(N,t2,neighbour_vec_s2,len_n_s2,-1,rand_cluster_degree,-1,-1);
	update_clustering_single_step(N,s1,neighbour_vec_t1,len_n_t1,-1,rand_cluster_degree,-1,-1);
	update_clustering_single_step(N,s2,neighbour_vec_t2,len_n_t2,-1,rand_cluster_degree,-1,-1);


	/* increment by 1/(k 2) rand_cluster_degree for nodes of classes 2,4,6,8
	   note that same node can be incremented and then decremented. don't change for s1--s2 or t1--t2 */
	update_clustering_single_step(N,t2,neighbour_vec_s1,len_n_s1,1,rand_cluster_degree,s2,t1);
	update_clustering_single_step(N,t1,neighbour_vec_s2,len_n_s2,1,rand_cluster_degree,s1,t2);
	update_clustering_single_step(N,s2,neighbour_vec_t1,len_n_t1,1,rand_cluster_degree,t2,s1);
	update_clustering_single_step(N,s1,neighbour_vec_t2,len_n_t2,1,rand_cluster_degree,t1,s2);

	if (neighbour_vec_s1!=NULL)
		free(neighbour_vec_s1);
	if (neighbour_vec_s2!=NULL)
		free(neighbour_vec_s2);
	if (neighbour_vec_t1!=NULL)
		free(neighbour_vec_t1);
	if (neighbour_vec_t2!=NULL)
		free(neighbour_vec_t2);

}

/* increments or decrements rand_cluster_degree according to sign for each neighbour which is connected to node
   unless the neighbour is node_exclude */
void update_clustering_single_step(Network *N,int node,int *neighbours,int len,int sign,double *rand_cluster_degree,int node_exclude1,int node_exclude2)
{
	register int i;
	double k;
	double temp;
	double inc;

	for (i=0;i<len;i++)
	{
		if ((is_edge(N,node,neighbours[i])||is_edge(N,neighbours[i],node))&&(neighbours[i]!=node_exclude1)&&(neighbours[i]!=node_exclude2))
		{
			k=(double)((N->doubledeg[neighbours[i]])+(N->indeg[neighbours[i]])+(N->outdeg[neighbours[i]]));
			temp=rand_cluster_degree[neighbours[i]];
			inc=(double)sign/(double)(k*(k-1)/2);
			inc=(double)(inc)*(0.5);	// this is because the cluster degree of this node is changed twice, once when
			//calling the function with s and once with t.
			rand_cluster_degree[neighbours[i]]=rand_cluster_degree[neighbours[i]]+inc;
			// in addition update the clustering coefficient of node
			k=(double)((N->doubledeg[node])+(N->indeg[node])+(N->outdeg[node]));
			inc=(double)sign/(double)(k*(k-1)/2);
			rand_cluster_degree[node]=rand_cluster_degree[node]+inc;
		}

	}

}

void copy_vec(double *target,double *source,int len)
{
	register int i;

	for (i=0;i<=len;i++)
		target[i]=source[i];

}

