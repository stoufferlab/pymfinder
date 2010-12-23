/************************************************************************
*
*  File name: output.c
*
*  Description: Output functions
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
#include "role.h"
#include "motif_ids.h"
#include "clustering.h"





/*************************** Global variables ****************************/


/******************************* Externs *********************************/


extern int DEBUG_LEVEL;

extern char *input_network_fname;
extern Network *G_N;
extern Network *G_WN;
extern Gnrl_st GNRL_ST;
//result table
extern Res_tbl RES_TBL, RES_TBL_SUB;
extern FILE *out_fp,*log_fp;
extern time_t start_time, end_time;


extern list64 *final_res, *final_res_all;
extern list64 *res_sub_motif;

//roles final results
extern Role_res *Roles_final_res;
extern role_st *Role_trans;


//Grassberger
extern int died_colonies;
extern int over_populated_colonies;



/******************************* Declerations ****************************/




/******************************** Functions  ******************************/

void
dump_network(FILE *fp,Network *N)
{
	int i,j;

	//printf("name :%s\t vertices: %d\t edges:%d\n",N->name,N->vertices_num,N->edges_num);
	if(N->edges_num>0) {
		for(i=1;i<=N->edges_num;i++)
			fprintf(fp,"%d %d\t", N->e_arr[i].s, N->e_arr[i].t);
		fprintf(fp,"\n");
	}
	if(N->e_arr_dbl>0) {
		for(i=1;i<=N->e_dbl_num;i++)
			fprintf(fp,"%d %d\t", N->e_arr_dbl[i].s, N->e_arr_dbl[i].t);
		fprintf(fp,"\n");
	}
	if(N->e_sin_num>0) {
		for(i=1;i<=N->e_sin_num;i++)
			fprintf(fp,"%d %d\t", N->e_arr_sin[i].s, N->e_arr_sin[i].t);
		fprintf(fp,"\n");
	}

	for (i=1; i<=N->vertices_num; i++) {
		fprintf(fp,"\n");
		for (j=1; j<=N->vertices_num; j++)
			fprintf(fp,"%d ",MatGet(N->mat,i,j));
	}
	fprintf(fp,"\n");
}

void
dump_motif_matrix(FILE *fp, int64 mtf_id)
{
	Matrix *m;

	m=init_matrix(GNRL_ST.mtf_sz);
	fill_mat_id(m,mtf_id);
	dump_matrix(fp,m,"");
	free_matrix(m);
	free(m);
}



void
dump_motifs_res(list64 *res, char *network_name)
{
	list64_item *l_res_id;

	//printf("motifs results for %s :\n", network_name);
	for(l_res_id=list64_get_next(res,NULL); l_res_id!=NULL;
		l_res_id=list64_get_next(res,l_res_id)) {
			//fprintf(GNRL_ST.out_fp,"id : %d  %d times.\n",l_res_id->val,*(int*)l_res_id->p);
		}
}

void
dump_motif_roles_count_by_vector(FILE *fp,Motif_res *mtf_res)
{
	int role_id;
	int i;
	int mtf_id;
	Role_res *role_res;

	mtf_id=(int)mtf_res->id;
	fprintf(fp,"\n\tRoles:  ");
	for(i=1;i<=3;i++) {
		role_id=Role_trans[mtf_id].roles[i];
		role_res = &Roles_final_res[role_id];
		if(Roles_final_res[role_id].role_id != role_id)
			printf("Error with role_id translation\n");
		if(i!=1)
			fprintf(fp,"\t\t");
		fprintf(fp,"[%d]: %d (%.1f+-%.1f)\n",
			i,
			(int)role_res->real_count,
			role_res->rand_mean,
			role_res->rand_std_dev);
	}
	fprintf(fp,"\n");
}


void
dump_motif_roles_count(FILE *fp, Motif_res *mtf_res)
{
	list_item *l_i;
	Member *members_arr;
	list *role_1,*role_2,*role_3;
	int i;

	list_init(&role_1);
	list_init(&role_2);
	list_init(&role_3);

	for(l_i=list_get_next(mtf_res->all_members,NULL); l_i!=NULL;
	l_i=list_get_next(mtf_res->all_members,l_i)) {
		members_arr=(Member*)(l_i->p);
		//pass through members array
		for(i=0;i<GNRL_ST.mtf_sz;i++) {
			//fill role lists according to node role :
			switch(members_arr[i].role) {
			case 1:
				if(list_get(role_1,members_arr[i].node) == NULL)
					list_insert(role_1,members_arr[i].node,NULL);
				break;
			case 2:
				if(list_get(role_2,members_arr[i].node) == NULL)
					list_insert(role_2,members_arr[i].node,NULL);
				break;
			case 3:
				if(list_get(role_3,members_arr[i].node) == NULL)
					list_insert(role_3,members_arr[i].node,NULL);
				break;
			default:
				break;
			}
		}
	}
	//only one role for that  id
	if(role_2->size == 0)
		fprintf(fp,"\tRoles: 1: %d\n",role_1->size);
	//only two roles
	else if (role_3->size == 0)
		fprintf(fp,"\tRoles appearances: 1: %d , 2: %d\n",role_1->size,role_2->size);
		//three roles
		else
			fprintf(fp,"\tRoles appearances: 1: %d , 2: %d , 3: %d\n",role_1->size,role_2->size,role_3->size);
	list_free_mem(role_1);
	list_free_mem(role_2);
	list_free_mem(role_3);
}


//dump roles results of all subgraphs
void
dump_roles_ratio_results(FILE *fp, list64* full_res)
{
	list64_item *l_mtf;
	Motif_res *mtf_res;
	int i;
	int mtf_id,role_id;
	Role_res *role_res;

	fprintf(fp,"\n\n\n\n\n\tRoles Power Results\n");
	fprintf(fp,"\t====================\n\n");

	fprintf(fp,"MOTIF\tROLE\tROLE\tREAL\tRAND\t\tPOWER\tPOWER\n");
	fprintf(fp,"ID\t\tID\tPOWER\tPOWER\t\tZSCORE\tPVAL\t\n");
	for(l_mtf=list64_get_next(full_res,NULL); l_mtf!=NULL;l_mtf=list64_get_next(full_res,l_mtf)) {
		mtf_res=(Motif_res*)l_mtf->p;
		mtf_id=(int)mtf_res->id;
		for(i=1;i<=3;i++) {
			role_id=Role_trans[mtf_id].roles[i];
			role_res = &Roles_final_res[role_id];
			//only for relevant roles (1,2,3 or 1,2 etc)
			if(Role_trans[mtf_id].roles[i] != 0) {
				if(i==1)
					fprintf(fp,"%d\t",mtf_id);
				else
					fprintf(fp,"\t");
				fprintf(fp,"%d\t%d\t%.2f\t%.2f+-%.2f\t%.2f\t%.3f\n",
					i,
					role_id,
					(double)role_res->real_ratio,
					(double)role_res->rand_ratio_mean,
					(double)role_res->rand_ratio_std_dev,
					(double)role_res->real_ratio_zscore,
					(double)role_res->real_ratio_pval);
			}
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
}

//dump roles results of all subgraphs
void
dump_roles_indeg_outdeg_stats(FILE *fp, list64* full_res)
{
	list64_item *l_mtf;
	Motif_res *mtf_res;
	int i;
	int mtf_id,role_id;
	Role_res *role_res;

	fprintf(fp,"\n\n\n\n\n\tRoles incoming outgoing degree Statistics\n");
	fprintf(fp,"\t====================================================\n\n");

	fprintf(fp,"MOTIF\tROLE\tROLE\tOUTDEG\t\tOUTDEG\tINDEG\t\tINDEG\n");
	fprintf(fp,"ID\t\tID\tSTATS\t\tMAX\tSTATS\t\tMAX\t\n");
	for(l_mtf=list64_get_next(full_res,NULL); l_mtf!=NULL;l_mtf=list64_get_next(full_res,l_mtf)) {
		mtf_res=(Motif_res*)l_mtf->p;
		mtf_id=(int)mtf_res->id;
		for(i=1;i<=3;i++) {
			role_id=Role_trans[mtf_id].roles[i];
			role_res = &Roles_final_res[role_id];
			//only for relevant roles (1,2,3 or 1,2 etc)
			if(Role_trans[mtf_id].roles[i] != 0) {
				if(i==1)
					fprintf(fp,"%d\t",mtf_id);
				else
					fprintf(fp,"\t");
				fprintf(fp,"%d\t%d\t%.1f+-%.1f\t%d\t%.1f+-%.1f\t%d\n",
					i,
					role_id,
					(double)role_res->outdeg_mean,
					(double)role_res->outdeg_std_dev,
					(int)role_res->outdeg_max,
					(double)role_res->indeg_mean,
					(double)role_res->indeg_std_dev,
					(int)role_res->indeg_max
					);
			}
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
}





//dump roles results of all subgraphs
void
dump_roles_results(FILE *fp, list64* full_res,int rnd_net_num)
{
	list64_item *l_mtf;
	Motif_res *mtf_res;
	int i;
	int mtf_id,role_id;
	Role_res *role_res;

	fprintf(fp,"\tSUMMARY ROLES RESULT\n");
	fprintf(fp,"\t====================\n\n");

	fprintf(fp,"\tNetwork name: %s\n",G_N->name);
	fprintf(fp,"\tMotif size searched is %d\n",GNRL_ST.mtf_sz);
	fprintf(fp,"\tNumber of random networks generated is: %d\n\n",rnd_net_num);

	fprintf(fp,"\t( Total num of different roles in all motifs size %d is : %d )\n\n",GNRL_ST.mtf_sz,30);


	fprintf(fp,"MOTIF\tROLE\tROLE\tNREAL\tNRAND\t\tNREAL\tNREAL\n");
	fprintf(fp,"ID\t\tID\t\tSTATS\t\tZSCORE\tPVAL\t\n");
	for(l_mtf=list64_get_next(full_res,NULL); l_mtf!=NULL;l_mtf=list64_get_next(full_res,l_mtf)) {
		mtf_res=(Motif_res*)l_mtf->p;
		mtf_id=(int)mtf_res->id;
		for(i=1;i<=3;i++) {
			role_id=Role_trans[mtf_id].roles[i];
			role_res = &Roles_final_res[role_id];
			//only for relevant roles (1,2,3 or 1,2 etc)
			if(Role_trans[mtf_id].roles[i] != 0) {
				if(i==1)
					fprintf(fp,"%d\t",mtf_id);
				else
					fprintf(fp,"\t");
				fprintf(fp,"%d\t%d\t%d\t%.1f+-%.1f\t%.2f\t%.3f\n",
					i,
					role_id,
					(int)role_res->real_count,
					(double)role_res->rand_mean,
					(double)role_res->rand_std_dev,
					(double)role_res->real_zscore,
					(double)role_res->real_pval);
			}
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n\n");

	dump_roles_ratio_results(fp,full_res);
	dump_roles_indeg_outdeg_stats(fp, full_res);
}







//dump members of all subgraphs
void
dump_members_lists(FILE *fp,list64 *res)
{
	list64_item *l_mtf;
	Motif_res *mtf_res;
	list_item *l_members;
	int i,j;
	Member* members_arr;
	double creal;
	int nreal;

	for(l_mtf=list64_get_next(res,NULL); l_mtf!=NULL;
		l_mtf=list64_get_next(res,l_mtf)) {
		mtf_res=(Motif_res*)l_mtf->p;
		if(mtf_res->all_members->size > 0) {
			if(GNRL_ST.run_prob_app) {
				creal=mtf_res->conc_real;
				fprintf(fp,"subgraph id = %lld\n", (uint64)mtf_res->id);
				fprintf(fp,"Creal : %.3f\n", creal);
			}
			else {
				nreal=(int)mtf_res->real_count;
				fprintf(fp,"subgraph id = %lld\n", (uint64)mtf_res->id);
				fprintf(fp,"Nreal : %d\n", nreal);
			}
			fprintf(fp,"=================================\n");
			//check if full or partial lists of memebers
			if(mtf_res->all_members->size==(int)mtf_res->real_count)
				fprintf(fp,"Full list of %d members:\n",mtf_res->all_members->size);
			else
				fprintf(fp,"Partial list of %d members:\n",mtf_res->all_members->size);
			//should be commented for ver1.1-web

			//Nadav - this is removed from the web application to keep compatibility for mfinder 1.1
#if 0
			//print roles
			//temp - get first member and extract roles from the members arr
			l_members=list_get_next(mtf_res->all_members,NULL);
			members_arr=(Member*)l_members->p;
			for(i=1;i<=GNRL_ST.mtf_sz;i++) {
				for(j=0;j<GNRL_ST.mtf_sz;j++) {
					if(members_arr[j].role==i)
						fprintf(fp,"(%d)\t",members_arr[j].role);
				}
			}
			fprintf(fp,"\n");
#endif

			fprintf(fp,"\n");
			//print members according to role order
			for(l_members=list_get_next(mtf_res->all_members,NULL); l_members!=NULL;
			l_members=list_get_next(mtf_res->all_members,l_members)) {
				members_arr=(Member*)l_members->p;
				//print members according to role order
				for(i=1;i<=GNRL_ST.mtf_sz;i++) {
					for(j=0;j<GNRL_ST.mtf_sz;j++) {
						if(members_arr[j].role==i)
							fprintf(fp,"%d\t",members_arr[j].node);
					}
				}
				fprintf(fp,"\n");
			}
			fprintf(fp,"\n");
		}
	}
}




void
gen_motif_ids_list(list64*res,list64**motif_ids_list_p, int rand_net_num, int with_dangling)
{
	list64_item *l_mtf;
	Motif_res *mtf_res;
	list64 *motif_ids_list;

	list64_init(&motif_ids_list);
	if(with_dangling==TRUE) {
		for(l_mtf=list64_get_next(res,NULL);l_mtf!=NULL;l_mtf=list64_get_next(res,l_mtf)) {
			mtf_res=(Motif_res*)l_mtf->p;
			if(GNRL_ST.run_prob_app==FALSE){
				//deterministic
				if(rand_net_num<1000){
					//use zscore
					if( ((mtf_res->real_zscore > GNRL_ST.zfactor_th)
						&& (mtf_res->real_count >=mtf_res->rand_mean*GNRL_ST.mfactor_th)
						&& (mtf_res->unique_appear>=GNRL_ST.unique_th))
						||
						((mtf_res->real_count>1) && (mtf_res->rand_mean==0)
						&&
						(GNRL_ST.rnd_net_num>0)&& (mtf_res->unique_appear>=GNRL_ST.unique_th))){

						list64_insert(motif_ids_list,(int)mtf_res->id,NULL);
					}
				}else{
					//use pval
					if( (mtf_res->real_pval < GNRL_ST.pval_th)
						&& (mtf_res->real_count >=mtf_res->rand_mean*GNRL_ST.mfactor_th)
						&& (mtf_res->unique_appear>=GNRL_ST.unique_th)) {
						list64_insert(motif_ids_list,(int)mtf_res->id,NULL);
					}
				}
			}else{
				//probabilisitc
				if(rand_net_num<1000){
					//use zscore

					if( ((mtf_res->conc_real_zscore > GNRL_ST.zfactor_th)
						&& (mtf_res->conc_real >=mtf_res->conc_rand_mean*GNRL_ST.mfactor_th)
						&& (mtf_res->unique_appear>=GNRL_ST.unique_th))
						||
						((mtf_res->conc_real>0) && (mtf_res->conc_rand_mean==0)
						&& (GNRL_ST.rnd_net_num>0) && (mtf_res->unique_appear>=GNRL_ST.unique_th)) ){
						list64_insert(motif_ids_list,(int)mtf_res->id,NULL);
					}
				}else{
					//use pval
					if( (mtf_res->conc_real_pval < GNRL_ST.pval_th)
						&& (mtf_res->conc_real >=mtf_res->conc_rand_mean*GNRL_ST.mfactor_th)
						&& (mtf_res->unique_appear>=GNRL_ST.unique_th)) {
						list64_insert(motif_ids_list,(int)mtf_res->id,NULL);
					}
				}
			}
		}
	}else{
		//Without dangling
		for(l_mtf=list64_get_next(res,NULL);l_mtf!=NULL;l_mtf=list64_get_next(res,l_mtf)) {
			mtf_res=(Motif_res*)l_mtf->p;
			if(mtf_res->dangling==TRUE)
				continue;
			if(GNRL_ST.run_prob_app==FALSE){
				//deterministic
				if(rand_net_num<1000){
					//use zscore
					if( ((mtf_res->real_zscore > GNRL_ST.zfactor_th)
						&& (mtf_res->real_count >=mtf_res->rand_mean*GNRL_ST.mfactor_th)
						&& (mtf_res->unique_appear>=GNRL_ST.unique_th))
						||
						((mtf_res->real_count>1) && (mtf_res->rand_mean==0)
						&&
						(GNRL_ST.rnd_net_num>0)&& (mtf_res->unique_appear>=GNRL_ST.unique_th))){

						list64_insert(motif_ids_list,(int)mtf_res->id,NULL);
					}
				}else{
					//use pval
					if( (mtf_res->real_pval < GNRL_ST.pval_th)
						&& (mtf_res->real_count >=mtf_res->rand_mean*GNRL_ST.mfactor_th)
						&& (mtf_res->unique_appear>=GNRL_ST.unique_th)) {
						list64_insert(motif_ids_list,(int)mtf_res->id,NULL);
					}
				}
			}else{
				//probabilisitc
				if(rand_net_num<1000){
					//use zscore

					if( ((mtf_res->conc_real_zscore > GNRL_ST.zfactor_th)
						&& (mtf_res->conc_real >=mtf_res->conc_rand_mean*GNRL_ST.mfactor_th)
						&& (mtf_res->unique_appear>=GNRL_ST.unique_th))
						||
						((mtf_res->conc_real>0) && (mtf_res->conc_rand_mean==0)
						&& (GNRL_ST.rnd_net_num>0) && (mtf_res->unique_appear>=GNRL_ST.unique_th)) ){
						list64_insert(motif_ids_list,(int)mtf_res->id,NULL);
					}
				}else{
					//use pval
					if( (mtf_res->conc_real_pval < GNRL_ST.pval_th)
						&& (mtf_res->conc_real >=mtf_res->conc_rand_mean*GNRL_ST.mfactor_th)
						&& (mtf_res->unique_appear>=GNRL_ST.unique_th)) {
						list64_insert(motif_ids_list,(int)mtf_res->id,NULL);
					}
				}
			}
		}
	}

	*motif_ids_list_p=motif_ids_list;
}

void
gen_top_zscore_list(list64*res,list64**top_zscore_list_p,list64 *motif_ids)
{
	list64_item *l_mtf,*l_top, *l_id;
	Motif_res *mtf_res;
	list64 *top_zscore_ids;
	double min_zscore=INFINITY;
	double *zscore;

	list64_init(&top_zscore_ids);

	for(l_id=list64_get_next(motif_ids,NULL); l_id!=NULL;l_id=list64_get_next(motif_ids,l_id)) {
		l_mtf=list64_get(res,l_id->val);
		mtf_res=(Motif_res*)l_mtf->p;
		if(GNRL_ST.run_prob_app==FALSE){
			//full enumeration
			if(top_zscore_ids->size<GNRL_ST.top_motifs_list_sz){
				if(mtf_res->dangling==FALSE) {
					zscore=(double*)calloc(1,sizeof(double));
					*zscore=mtf_res->real_zscore;
					list64_insert(top_zscore_ids,(int)mtf_res->id,(void*)zscore);
					//zscore of appearances
					if(mtf_res->real_zscore<min_zscore)
						min_zscore=mtf_res->real_zscore;
				}
			}else{ //bigger then TOP_MOTIF_LIST_SZ , delete the min weight and insert the new one
				if((mtf_res->real_zscore > min_zscore) && (mtf_res->dangling==FALSE)){
					for(l_top=list64_get_next(top_zscore_ids,NULL);
					l_top!=NULL;l_top=list64_get_next(top_zscore_ids,l_top)) {
						if(*(double*)(l_top->p)==min_zscore){
							list64_delete(top_zscore_ids,l_top->val);
							break;
						}

					}
					zscore=(double*)calloc(1,sizeof(double));
					*zscore=mtf_res->real_zscore;
					list64_insert(top_zscore_ids,(int)mtf_res->id,(void*)zscore);
					//update min_zscore
					min_zscore=INFINITY;
					for(l_top=list64_get_next(top_zscore_ids,NULL);
					l_top!=NULL;l_top=list64_get_next(top_zscore_ids,l_top)) {
						if(*(double*)(l_top->p)<min_zscore)
							min_zscore=*(double*)(l_top->p);
					}
				}
			}
		}else{
			//probabilistic alg
			if(top_zscore_ids->size<GNRL_ST.top_motifs_list_sz){
				if(mtf_res->dangling==FALSE) {
					zscore=(double*)calloc(1,sizeof(double));
					*zscore=mtf_res->conc_real_zscore;
					list64_insert(top_zscore_ids,(int)mtf_res->id,(void*)zscore);
					//zscore of apperance
					if(mtf_res->conc_real_zscore<min_zscore)
						min_zscore=mtf_res->conc_real_zscore;
				}
			}else{ //bigger then GNRL_ST.top_motifs_list_sz , delete thge min weight and insert the new one
				if((mtf_res->conc_real_zscore > min_zscore)&& (mtf_res->dangling==FALSE)){
					for(l_top=list64_get_next(top_zscore_ids,NULL);
					l_top!=NULL;l_top=list64_get_next(top_zscore_ids,l_top)) {
						if(*(double*)(l_top->p)==min_zscore){
							list64_delete(top_zscore_ids,l_top->val);
							break;
						}
					}
					zscore=(double*)calloc(1,sizeof(double));
					*zscore=mtf_res->conc_real_zscore;
					list64_insert(top_zscore_ids,(int)mtf_res->id,(void*)zscore);
					//update min_zscore
					min_zscore=INFINITY;
					for(l_top=list64_get_next(top_zscore_ids,NULL);
					l_top!=NULL;l_top=list64_get_next(top_zscore_ids,l_top)) {
						if(*(double*)(l_top->p)<min_zscore)
							min_zscore=*(double*)(l_top->p);
					}
				}
			}
		}
	}
	*top_zscore_list_p=top_zscore_ids;
}


void
gen_top_conc_list(list64*res,list64**top_conc_list_p,list64 *motif_ids)
{
	list64_item *l_mtf,*l_top, *l_id;
	Motif_res *mtf_res;
	list64 *top_conc_ids;
	double min_conc=INFINITY;
	double *conc;

	list64_init(&top_conc_ids);

	for(l_id=list64_get_next(motif_ids,NULL); l_id!=NULL;l_id=list64_get_next(motif_ids,l_id)) {
		l_mtf=list64_get(res,l_id->val);
		mtf_res=(Motif_res*)l_mtf->p;
		if(GNRL_ST.run_prob_app==FALSE){
			//full enumeration
			if(top_conc_ids->size<GNRL_ST.top_motifs_list_sz){
				//if not dangling then insert
				if(mtf_res->dangling==FALSE) {
					conc=(double*)calloc(1,sizeof(double));
					*conc=mtf_res->conc_real;
					list64_insert(top_conc_ids,(int)mtf_res->id,(void*)conc);
					//zscore of appearances
					if(mtf_res->conc_real<min_conc)
						min_conc=mtf_res->conc_real;
				}
			}else{ //bigger then GNRL_ST.top_motifs_list_sz , delete the min weight and insert the new one
				if( (mtf_res->conc_real > min_conc) && (mtf_res->dangling==FALSE) ){
					for(l_top=list64_get_next(top_conc_ids,NULL);
					l_top!=NULL;l_top=list64_get_next(top_conc_ids,l_top)) {
						if(*(double*)(l_top->p)==min_conc){
							list64_delete(top_conc_ids,l_top->val);
							break;
						}

					}
					conc=(double*)calloc(1,sizeof(double));
					*conc=mtf_res->conc_real;
					list64_insert(top_conc_ids,(int)mtf_res->id,(void*)conc);
					//update min_conc
					min_conc=INFINITY;
					for(l_top=list64_get_next(top_conc_ids,NULL);
					l_top!=NULL;l_top=list64_get_next(top_conc_ids,l_top)) {
						if(*(double*)(l_top->p)<min_conc)
							min_conc=*(double*)(l_top->p);
					}
				}
			}
		}else{
			//probabilistic alg
			if(top_conc_ids->size<GNRL_ST.top_motifs_list_sz){
				//if not dangling then insert
				if(mtf_res->dangling==FALSE) {
					conc=(double*)calloc(1,sizeof(double));
					*conc=mtf_res->conc_real;
					list64_insert(top_conc_ids,(int)mtf_res->id,(void*)conc);
					//conc of apperance
					if(mtf_res->conc_real<min_conc)
						min_conc=mtf_res->conc_real;
				}
			}else{ //bigger then TOP_MOTIF_LIST_SZ , delete thge min weight and insert the new one
				if( (mtf_res->conc_real > min_conc) && (mtf_res->dangling==FALSE) ){
					for(l_top=list64_get_next(top_conc_ids,NULL);
					l_top!=NULL;l_top=list64_get_next(top_conc_ids,l_top)) {
						if(*(double*)(l_top->p)==min_conc){
							list64_delete(top_conc_ids,l_top->val);
							break;
						}
					}
					conc=(double*)calloc(1,sizeof(double));
					*conc=mtf_res->conc_real;
					list64_insert(top_conc_ids,(int)mtf_res->id,(void*)conc);
					//update min_conc
					min_conc=INFINITY;
					for(l_top=list64_get_next(top_conc_ids,NULL);
					l_top!=NULL;l_top=list64_get_next(top_conc_ids,l_top)) {
						if(*(double*)(l_top->p)<min_conc)
							min_conc=*(double*)(l_top->p);
					}
				}
			}
		}
	}
	*top_conc_list_p=top_conc_ids;
}

void
gen_top_motifs_list(list64*res,list64**top_motifs_list_p,list64 *motif_ids)
{
	list64_item *l_mtf,*l_top, *l_id;
	Motif_res *mtf_res;
	list64 *top_ids;
	double min_score=INFINITY;
	double *score;
	double score_val;

	list64_init(&top_ids);

	for(l_id=list64_get_next(motif_ids,NULL); l_id!=NULL;l_id=list64_get_next(motif_ids,l_id)) {
		l_mtf=list64_get(res,l_id->val);
		mtf_res=(Motif_res*)l_mtf->p;
		//score=concentration *2*log10(zscore)
		if(GNRL_ST.run_prob_app==FALSE)
			score_val=mtf_res->conc_real*(2*log10(mtf_res->real_zscore));
		else
			score_val=mtf_res->conc_real*(2*log10(mtf_res->conc_real_zscore));

		if(top_ids->size<GNRL_ST.top_motifs_list_sz){
			score=(double*)calloc(1,sizeof(double));
			*score=score_val;
			list64_insert(top_ids,(int)mtf_res->id,(void*)score);
			//conc of apperance
			if(score_val<min_score)
				min_score=score_val;
		}else{ //bigger then GNRL_ST.top_motifs_list_sz , delete thge min weight and insert the new one
			if(score_val > min_score){
				for(l_top=list64_get_next(top_ids,NULL);
				l_top!=NULL;l_top=list64_get_next(top_ids,l_top)) {
					if(*(double*)(l_top->p)==min_score){
						list64_delete(top_ids,l_top->val);
						break;
					}
				}
				score=(double*)calloc(1,sizeof(double));
				*score=score_val;
				list64_insert(top_ids,(int)mtf_res->id,(void*)score);
				//update min_conc
				min_score=INFINITY;
				for(l_top=list64_get_next(top_ids,NULL);
				l_top!=NULL;l_top=list64_get_next(top_ids,l_top)) {
					if(*(double*)(l_top->p)<min_score)
						min_score=*(double*)(l_top->p);
				}
			}
		}
	}
	*top_motifs_list_p=top_ids;
}



void
dump_network_stat(FILE *fp,Network *N, int rnd_net_num, list64 *res)
{
	list64_item *l_i;
	Motif_res *mtf_res;
	int subgraphs_total_num=0;

	fprintf(fp,"\nMOTIF FINDER RESULTS:\n\n");
	fprintf(fp,"\tNetwork name: %s\n",G_N->name);
	if(GNRL_ST.undirected_flag==FALSE) {
		fprintf(fp,"\tNetwork type: Directed\n");
		fprintf(fp,"\tNum of Nodes: %d Num of Edges: %d\n", G_N->vertices_num, G_N->edges_num);
	}else{
		fprintf(fp,"\tNetwork type: Non-Directed\n");
		fprintf(fp,"\tNum of Nodes: %d Num of Edges: %d\n", G_N->vertices_num, G_N->edges_num/2);
	}
	fprintf(fp,"\tNum of Nodes with edges: %d\n", G_N->con_vertices_num);
	if(GNRL_ST.calc_self_edges == TRUE)
	{
		int num_of_self_edges = 0, i;

		for(i = 0; i < G_N->vertices_num; i++)
		{
			num_of_self_edges += G_N->self_edge[i];
		}
		fprintf(fp,"\tNum of Nodes with self edges: %d\n", num_of_self_edges);
	}
	//calc total num of subgraphs
	if( (GNRL_ST.run_prob_app==FALSE) && (res!=NULL) ) {
		for(l_i=list64_get_next(res,NULL);l_i!=NULL;l_i=list64_get_next(res,l_i)){
			mtf_res=(Motif_res*)l_i->p;
			subgraphs_total_num+=(int)mtf_res->real_count;
		}
	}

	if((GNRL_ST.quiet_mode==FALSE)) {

		fprintf(fp,"\tMaximal out degree (out-hub) : %d\n",G_N->out_hub_deg);
		fprintf(fp,"\tMaximal in degree (in-hub) : %d\n",G_N->in_hub_deg);
		fprintf(fp,"\tRoots num: %d Leaves num: %d\n", G_N->roots_num, G_N->leafs_num);
		fprintf(fp,"\tSingle Edges num: %d Mutual Edges num: %d\n", G_N->e_sin_num,G_N->e_dbl_num/2);
		fprintf(fp,"\n\tMotif size searched  %d\n",GNRL_ST.mtf_sz);
		if( (GNRL_ST.run_prob_app==FALSE) && (subgraphs_total_num!=0) )
			fprintf(fp,"\tTotal number of %d-node subgraphs : %d\n",GNRL_ST.mtf_sz,subgraphs_total_num);
		fprintf(fp,"\tNumber of random networks generated : %d\n",rnd_net_num);
		if(GNRL_ST.use_clustering)
			fprintf(fp,"\tRandom networks generation method: Preserve clustering sequence\n");
		else if(GNRL_ST.use_stubs_method)
			fprintf(fp,"\tRandom networks generation method: Stubs\n");
		else if(GNRL_ST.use_metropolis)
			fprintf(fp,"\tRandom networks generation method: Switches (with metropolis)\n");
		else if(GNRL_ST.r_grassberger){
			fprintf(fp,"\tRandom networks generation method: Grassberger\n");
			fprintf(fp,"\tColony size: %d, Reff: %.2f Clones: %.3f+-%.3f\n",
				GNRL_ST.r_grass_colony_sz,GNRL_ST.rnstat.grass_reff,
				GNRL_ST.rnstat.clones_mean,GNRL_ST.rnstat.clones_stddev);
			fprintf(fp,"\tDied colonies :%d Over populated colonies %d\n",died_colonies,over_populated_colonies);
		}else{
			fprintf(fp,"\tRandom networks generation method: Switches\n");
			fprintf(fp,"\tNum of Switches range: %.1f-%.1f, Success switches Ratio:%.3f+-%.2f\n",
				GNRL_ST.r_switch_factor,(double)(2*GNRL_ST.r_switch_factor),
				GNRL_ST.rnstat.switch_ratio_mean,GNRL_ST.rnstat.switch_ratio_stddev);
		}
	}
	fprintf(fp,"\n");
}






//dump results to screen and to file pointer (fp parameter)
void
dump_final_res(FILE *fp,list64 *res, list64* full_res, int rnd_net_num)
{
	list64_item *l_mtf,*l_id;
	Motif_res *mtf_res;
	list64 *motif_ids; //list of all subgraphs passed criteria for motifs
	list64 *non_dangling_motif_ids; //list of all non-dangling subgraphs passed criteria for motifs
	list64 *top_motif_ids; //list of top motifs by concentration

	//Generate motif ids list
	//according to motifs test criteria

	//with dangling edges
	gen_motif_ids_list(res,&motif_ids, rnd_net_num, 1);
	//without dangling edges
	gen_motif_ids_list(res,&non_dangling_motif_ids, rnd_net_num, 0);

	//if more then MAX_FULL_MOTIF_LIST_SZ (default=10) motifs were found then
	//make lists of top in Zscore and top in Conc
	if( motif_ids->size>=MAX_FULL_MOTIF_LIST_SZ) {

		gen_top_motifs_list(res,&top_motif_ids,non_dangling_motif_ids);
	}

	if((GNRL_ST.quiet_mode==FALSE)) {

		// to stdout
		if(GNRL_ST.run_prob_app==FALSE) {
			//deterministic alg
			dump_network_stat(stdout,G_N,rnd_net_num,res);

			fprintf(stdout,"The following motifs were found:\n\n");
			//if num of networks >=1000 then use pval otherwise ignore it in the motifs criteria
			if(GNRL_ST.rnd_net_num<1000) {
				fprintf(stdout,"Criteria taken : Nreal Zscore > %.2f\n",GNRL_ST.zfactor_th);
				fprintf(stdout,"                 Pval ignored (due to small number of random networks)\n");
				fprintf(stdout,"                 Mfactor > %.2f\n",GNRL_ST.mfactor_th);
				if(GNRL_ST.calc_unique_flag)
					fprintf(stdout,"                 Uniqueness >= %d\n\n",GNRL_ST.unique_th);
				else
					fprintf(stdout,"                 Uniquness ignored\n\n");

				//dump top motif lists
				if(motif_ids->size>MAX_FULL_MOTIF_LIST_SZ) {
					fprintf(stdout,"\n\tTop Motifs List\n\n");
					fprintf(stdout,"MOTIF\tNREAL\tNRAND\t\tNREAL\tNREAL\tUNIQ\tCREAL\n");
					fprintf(stdout,"ID\t\tSTATS\t\tZSCORE\tPVAL\tVAL\t[MILI]\t\n\n");
					for(l_id=list64_get_next(top_motif_ids,NULL); l_id!=NULL;
					l_id=list64_get_next(top_motif_ids,l_id)) {
						l_mtf=list64_get(res,l_id->val);
						mtf_res=(Motif_res*)l_mtf->p;

						if(mtf_res->real_zscore!=UNDEFINED)
						fprintf(stdout,"%lld\t%.0f\t%.1f+-%.1f\t%.2f\t%.3f\t%d\t%.2f\n",
							(uint64)mtf_res->id,
							mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
							mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
							mtf_res->conc_real);
						else
							fprintf(stdout,"%lld\t%.0f\t%.1f+-%.1f\t%.0f\t%.3f\t%d\t%.2f\n",
							(uint64)mtf_res->id,
							mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
							mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
							mtf_res->conc_real);


						dump_motif_matrix(stdout,mtf_res->id);
						if(GNRL_ST.calc_roles==TRUE)
							dump_motif_roles_count_by_vector(stdout,mtf_res);
						fprintf(stdout,"\n");
					}
					fprintf(stdout,"\n\n\tTotal No. of non-dangling motifs : %d  \n\n",non_dangling_motif_ids->size);
					fprintf(stdout,"\n\n\tFull list includes %d motifs (listed in the output file) \n",motif_ids->size);
				}else{
					//less then MAX_FULL_MOTIF_LIST_SZ motifs output all the list
					fprintf(stdout,"\n\n\tFull list includes %d motifs\n\n",motif_ids->size);
					fprintf(stdout,"MOTIF\tNREAL\tNRAND\t\tNREAL\tNREAL\tUNIQ\tCREAL\n");
					fprintf(stdout,"ID\t\tSTATS\t\tZSCORE\tPVAL\tVAL\t[MILI]\t\n\n");
					for(l_id=list64_get_next(motif_ids,NULL); l_id!=NULL;
					l_id=list64_get_next(motif_ids,l_id)) {
						l_mtf=list64_get(res,l_id->val);
						mtf_res=(Motif_res*)l_mtf->p;
						if(mtf_res->real_zscore!=UNDEFINED)
						fprintf(stdout,"%lld\t%.0f\t%.1f+-%.1f\t%.2f\t%.3f\t%d\t%.2f\n",
							(uint64)mtf_res->id,
							mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
							mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
							mtf_res->conc_real);
						else
							fprintf(stdout,"%lld\t%.0f\t%.1f+-%.1f\t%.0f\t%.3f\t%d\t%.2f\n",
							(uint64)mtf_res->id,
							mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
							mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
							mtf_res->conc_real);
						dump_motif_matrix(stdout,mtf_res->id);
						if(GNRL_ST.calc_roles==TRUE)
							dump_motif_roles_count_by_vector(stdout,mtf_res);
						fprintf(stdout,"\n");
					}
				}
				//enough random network num - consider pval
			}else{
				fprintf(stdout,"Criteria taken : Pval < %.3f\n",GNRL_ST.pval_th);
				fprintf(stdout,"                 Mfactor > %.2f\n",GNRL_ST.mfactor_th);
				if(GNRL_ST.calc_unique_flag)
					fprintf(stdout,"                 Uniqueness >= %d\n\n",GNRL_ST.unique_th);
				else
					fprintf(stdout,"                 Uniquness ignored\n\n");

				fprintf(stdout,"MOTIF\tNREAL\tNRAND\t\tNREAL\tNREAL\tUNIQ\tCREAL\n");
				fprintf(stdout,"ID\t\tSTATS\t\tZSCORE\tPVAL\tVAL\t[MILI]\t\n");

				for(l_id=list64_get_next(motif_ids,NULL); l_id!=NULL;
				l_id=list64_get_next(motif_ids,l_id)) {
					l_mtf=list64_get(res,l_id->val);
					mtf_res=(Motif_res*)l_mtf->p;
					if(mtf_res->real_zscore!=UNDEFINED)
					fprintf(stdout,"%lld\t%.0f\t%.1f+-%.1f\t%.2f\t%.3f\t%d\t%.2f\n",
						(uint64)mtf_res->id,
						mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
						mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
						mtf_res->conc_real);
					else
						fprintf(stdout,"%lld\t%.0f\t%.1f+-%.1f\t%.0f\t%.3f\t%d\t%.2f\n",
						(uint64)mtf_res->id,
						mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
						mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
						mtf_res->conc_real);
					dump_motif_matrix(stdout,mtf_res->id);
					if(GNRL_ST.calc_roles==TRUE)
						dump_motif_roles_count_by_vector(stdout,mtf_res);
					fprintf(stdout,"\n");
				}
			}
			fprintf(stdout,"\n");

		} else {
			//probabilistic
			dump_network_stat(stdout,G_N,rnd_net_num,res);

			if(GNRL_ST.prob_converge_mode)
				fprintf(stdout,"\tAlgorithm used : Probabilistic-Convergness\n\t(real network num of samples : %d)\n\n",
				GNRL_ST.prob_total_samples_num[0]);
			else
				fprintf(stdout,"\tAlgorithm used : Probabilistic (num of samples : %d)\n\n",GNRL_ST.prob_base_samples_num);
			fprintf(stdout,"The following motifs were found:\n\n");


			//dump motifs according to concentrations

			//if num of networks >=1000 then use pval otherwise ignore it in the motifs criteria
			if(GNRL_ST.rnd_net_num<1000) {
				fprintf(stdout,"Criteria taken : Creal Zscore > %.2f\n",GNRL_ST.zfactor_th);
				fprintf(stdout,"                 Pval ignored (due to small number of random networks)\n");
				fprintf(stdout,"                 Mfactor > %.2f\n",GNRL_ST.mfactor_th);
				if(GNRL_ST.calc_unique_flag)
					fprintf(stdout,"                 Uniqueness >= %d\n\n",GNRL_ST.unique_th);
				else
					fprintf(stdout,"                 Uniquness ignored\n\n");


				//dump top motif lists
				if(motif_ids->size>MAX_FULL_MOTIF_LIST_SZ) {
					fprintf(stdout,"\n\tTop Motifs List\n\n");
					fprintf(stdout,"MOTIF\tCREAL\tCRAND\t\tCREAL\tCREAL\tUNIQ\tHITS\n");
					fprintf(stdout,"ID\t[mili]\tSTATS\t\tZSCORE\tPVAL\t\n\n");
					for(l_id=list64_get_next(top_motif_ids,NULL); l_id!=NULL;
					l_id=list64_get_next(top_motif_ids,l_id)) {
						l_mtf=list64_get(res,l_id->val);
						mtf_res=(Motif_res*)l_mtf->p;
						if(mtf_res->conc_real_zscore!=UNDEFINED)
						fprintf(stdout,"%lld\t%.3f\t%.3f+-%.3f\t%.2f\t%.3f\t%d\t%d\t\n",
							(uint64)mtf_res->id,
							mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
							mtf_res->conc_real_zscore,mtf_res->conc_real_pval,
							(int)mtf_res->unique_appear,(int)mtf_res->hits_num);
						else
							fprintf(stdout,"%lld\t%.3f\t%.3f+-%.3f\t%.0f\t%.3f\t%d\t%d\t\n",
							(uint64)mtf_res->id,
							mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
							mtf_res->conc_real_zscore,mtf_res->conc_real_pval,
							(int)mtf_res->unique_appear,(int)mtf_res->hits_num);
						dump_motif_matrix(stdout,mtf_res->id);
						fprintf(stdout,"\n");
					}
					fprintf(stdout,"\n\n\tTotal No. of non-dangling motifs : %d  \n\n",non_dangling_motif_ids->size);
					fprintf(stdout,"\n\n\tFull list includes %d motifs (listed in the output file) \n",motif_ids->size);
				}else{
					//less then MAX_FULL_MOTIF_LIST_SZ motifs output all the list
					fprintf(stdout,"\n\n\tFull list includes %d motifs\n\n",motif_ids->size);
					fprintf(stdout,"MOTIF\tCREAL\tCRAND\t\tCREAL\tCREAL\tUNIQ\tHITS\n");
					fprintf(stdout,"ID\t[mili]\tSTATS\t\tZSCORE\tPVAL\t\n\n");
					for(l_id=list64_get_next(motif_ids,NULL); l_id!=NULL;
					l_id=list64_get_next(motif_ids,l_id)) {
						l_mtf=list64_get(res,l_id->val);
						mtf_res=(Motif_res*)l_mtf->p;
						if(mtf_res->conc_real_zscore!=UNDEFINED)
						fprintf(stdout,"%lld\t%.3f\t%.3f+-%.3f\t%.2f\t%.3f\t%d\t%d\t\n",
							(uint64)mtf_res->id,
							mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
							mtf_res->conc_real_zscore,mtf_res->conc_real_pval,
							(int)mtf_res->unique_appear,(int)mtf_res->hits_num);
						else
							fprintf(stdout,"%lld\t%.3f\t%.3f+-%.3f\t%.0f\t%.3f\t%d\t%d\t\n",
							(uint64)mtf_res->id,
							mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
							mtf_res->conc_real_zscore,mtf_res->conc_real_pval,
							(int)mtf_res->unique_appear,(int)mtf_res->hits_num);
						dump_motif_matrix(stdout,mtf_res->id);
						fprintf(stdout,"\n");
					}
				}
			} else {
				//num of rand networks>=1000
				//consider pval
				fprintf(stdout,"Criteria taken : Pval < %.3f\n",GNRL_ST.pval_th);
				fprintf(stdout,"                 Mfactor > %.2f\n",GNRL_ST.mfactor_th);
				if(GNRL_ST.calc_unique_flag)
					fprintf(stdout,"                 Uniqueness >= %d\n\n",GNRL_ST.unique_th);
				else
					fprintf(stdout,"                 Uniquness ignored\n\n");

				fprintf(stdout,"MOTIF\tCREAL\tCRAND\t\tCREAL\tCREAL\tUNIQ\tHITS\n");
				fprintf(stdout,"ID\t[mili]\tSTATS\t\tZSCORE\tPVAL\t\n");

				for(l_id=list64_get_next(motif_ids,NULL); l_id!=NULL;
				l_id=list64_get_next(motif_ids,l_id)) {
					l_mtf=list64_get(res,l_id->val);
					mtf_res=(Motif_res*)l_mtf->p;
					if(mtf_res->conc_real_zscore!=UNDEFINED)
					fprintf(stdout,"%lld\t%.3f\t%.3f+-%.3f\t%.2f\t%.3f\t%d\t%d\t\n",
						(uint64)mtf_res->id,
						mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
						mtf_res->conc_real_zscore,mtf_res->conc_real_pval,
						(int)mtf_res->unique_appear,(int)mtf_res->hits_num);
					else
						fprintf(stdout,"%lld\t%.3f\t%.3f+-%.3f\t%.0f\t%.3f\t%d\t%d\t\n",
						(uint64)mtf_res->id,
						mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
						mtf_res->conc_real_zscore,mtf_res->conc_real_pval,
						(int)mtf_res->unique_appear,(int)mtf_res->hits_num);

					dump_motif_matrix(stdout,mtf_res->id);
					fprintf(stdout,"\n");
				}
			}
			fprintf(stdout,"\n");
		}
	}

	//dump motif results  to result file

	fprintf(fp,"mfinder Version %.2f\n",VERSION);

	if(GNRL_ST.run_prob_app==FALSE) {
		//deterministic alg
		dump_network_stat(fp,G_N,rnd_net_num,res);

		fprintf(fp,"The following motifs were found:\n\n");
		//if num of networks >=1000 then use pval otherwise ignore it in the motifs criteria
		if(GNRL_ST.rnd_net_num<1000) {
			fprintf(fp,"Criteria taken : Nreal Zscore > %.2f\n",GNRL_ST.zfactor_th);
			fprintf(fp,"                 Pval ignored (due to small number of random networks)\n");
			fprintf(fp,"                 Mfactor > %.2f\n",GNRL_ST.mfactor_th);
			if(GNRL_ST.calc_unique_flag)
				fprintf(fp,"                 Uniqueness >= %d\n\n",GNRL_ST.unique_th);
			else
				fprintf(fp,"                 Uniquness ignored\n\n");
			if(motif_ids->size>MAX_FULL_MOTIF_LIST_SZ) {
				fprintf(fp,"\n\tTop Motifs List\n\n");
				fprintf(fp,"MOTIF\tNREAL\tNRAND\t\tNREAL\tNREAL\tUNIQ\tCREAL\n");
				fprintf(fp,"ID\t\tSTATS\t\tZSCORE\tPVAL\tVAL\t[MILI]\t\n\n");
				for(l_id=list64_get_next(top_motif_ids,NULL); l_id!=NULL;
				l_id=list64_get_next(top_motif_ids,l_id)) {
					l_mtf=list64_get(res,l_id->val);
					mtf_res=(Motif_res*)l_mtf->p;
					if(mtf_res->real_zscore!=UNDEFINED)
					fprintf(fp,"%lld\t%.0f\t%.1f+-%.1f\t%.2f\t%.3f\t%d\t%.2f\n",
						(uint64)mtf_res->id,
						mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
						mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
						mtf_res->conc_real);
					else
						fprintf(fp,"%lld\t%.0f\t%.1f+-%.1f\t%.0f\t%.3f\t%d\t%.2f\n",
						(uint64)mtf_res->id,
						mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
						mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
						mtf_res->conc_real);
					dump_motif_matrix(fp,mtf_res->id);
					if(GNRL_ST.calc_roles==TRUE)
						dump_motif_roles_count_by_vector(stdout,mtf_res);
					fprintf(fp,"\n");
				}

				if(GNRL_ST.out_non_dangling_motifs) {
					fprintf(fp,"\n\n\tTotal No. of non-dangling motifs : %d\n\n",non_dangling_motif_ids->size);
					fprintf(fp,"MOTIF\tNREAL\tNRAND\t\tNREAL\tNREAL\tUNIQ\tCREAL\n");
					fprintf(fp,"ID\t\tSTATS\t\tZSCORE\tPVAL\tVAL\t[MILI]\t\n\n");
					for(l_id=list64_get_next(non_dangling_motif_ids,NULL); l_id!=NULL;
					l_id=list64_get_next(non_dangling_motif_ids,l_id)) {
						l_mtf=list64_get(res,l_id->val);
						mtf_res=(Motif_res*)l_mtf->p;
						if(mtf_res->real_zscore!=UNDEFINED)
						fprintf(fp,"%lld\t%.0f\t%.1f+-%.1f\t%.2f\t%.3f\t%d\t%.2f\n",
							(uint64)mtf_res->id,
							mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
							mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
							mtf_res->conc_real);
						else
							fprintf(fp,"%lld\t%.0f\t%.1f+-%.1f\t%.0f\t%.3f\t%d\t%.2f\n",
							(uint64)mtf_res->id,
							mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
							mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
							mtf_res->conc_real);
						dump_motif_matrix(fp,mtf_res->id);
						if(GNRL_ST.calc_roles==TRUE)
							dump_motif_roles_count_by_vector(stdout,mtf_res);
						fprintf(fp,"\n");
					}
				} else {
					fprintf(fp,"\n\n\t(Total No. of non-dangling motifs : %d)\n\n",non_dangling_motif_ids->size);
				}
			}
			fprintf(fp,"\n\n\tFull list includes %d motifs\n",motif_ids->size);
			fprintf(fp,"MOTIF\tNREAL\tNRAND\t\tNREAL\tNREAL\tUNIQ\tCREAL\n");
			fprintf(fp,"ID\t\tSTATS\t\tZSCORE\tPVAL\tVAL\t[MILI]\t\n\n");
			//dump all motifs
			for(l_id=list64_get_next(motif_ids,NULL); l_id!=NULL;
			l_id=list64_get_next(motif_ids,l_id)) {
				l_mtf=list64_get(res,l_id->val);
				mtf_res=(Motif_res*)l_mtf->p;

				if(mtf_res->real_zscore!=UNDEFINED)
				fprintf(fp,"%lld\t%.0f\t%.1f+-%.1f\t%.2f\t%.3f\t%d\t%.2f\n",
					(uint64)mtf_res->id,
					mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
					mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
					mtf_res->conc_real);
				else
					fprintf(fp,"%lld\t%.0f\t%.1f+-%.1f\t%.0f\t%.3f\t%d\t%.2f\n",
					(uint64)mtf_res->id,
					mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
					mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
					mtf_res->conc_real);
				dump_motif_matrix(fp,mtf_res->id);
				if(GNRL_ST.calc_roles==TRUE)
					dump_motif_roles_count_by_vector(fp,mtf_res);
				fprintf(fp,"\n");
			}
			//enough random network num - consider pval
		}else{
			fprintf(fp,"Criteria taken : Pval < %.3f\n",GNRL_ST.pval_th);
			fprintf(fp,"                 Mfactor > %.2f\n",GNRL_ST.mfactor_th);
			if(GNRL_ST.calc_unique_flag)
				fprintf(fp,"                 Uniqueness >= %d\n\n",GNRL_ST.unique_th);
			else
				fprintf(fp,"                 Uniquness ignored\n\n");

			fprintf(fp,"MOTIF\tNREAL\tNRAND\t\tNREAL\tNREAL\tUNIQ\tCREAL\n");
			fprintf(fp,"ID\t\tSTATS\t\tZSCORE\tPVAL\tVAL\t[MILI]\t\n");

			for(l_id=list64_get_next(motif_ids,NULL); l_id!=NULL;
			l_id=list64_get_next(motif_ids,l_id)) {
				l_mtf=list64_get(res,l_id->val);
				mtf_res=(Motif_res*)l_mtf->p;
				if(mtf_res->real_zscore!=UNDEFINED)
				fprintf(fp,"%lld\t%.0f\t%.1f+-%.1f\t%.2f\t%.3f\t%d\t%.2f\n",
					(uint64)mtf_res->id,
					mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
					mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
					mtf_res->conc_real);
				else
					fprintf(fp,"%lld\t%.0f\t%.1f+-%.1f\t%.0f\t%.3f\t%d\t%.2f\n",
					(uint64)mtf_res->id,
					mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
					mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
					mtf_res->conc_real);
				dump_motif_matrix(fp,mtf_res->id);
				if(GNRL_ST.calc_roles==TRUE)
					dump_motif_roles_count_by_vector(fp,mtf_res);
				fprintf(fp,"\n");
			}
		}
		fprintf(fp,"\n");

	} else {
		//probabilistic
		dump_network_stat(fp,G_N,rnd_net_num,res);
		if(GNRL_ST.prob_converge_mode)
			fprintf(fp,"\tAlgorithm used : Probabilistic-Convergness\n\n\t(real network num of samples : %d)\n\n",
			GNRL_ST.prob_total_samples_num[0]);
		else
			fprintf(fp,"\tAlgorithm used : Probabilistic (num of samples : %d)\n\n",GNRL_ST.prob_base_samples_num);
		fprintf(fp,"The following motifs were found:\n\n");


		//dump motifs according to concentrations

		//if num of networks >=1000 then use pval otherwise ignore it in the motifs criteria
		if(GNRL_ST.rnd_net_num<1000) {
			fprintf(fp,"Criteria taken : Creal Zscore > %.2f\n",GNRL_ST.zfactor_th);
			fprintf(fp,"                 Pval ignored (due to small number of random networks)\n");
			fprintf(fp,"                 Mfactor > %.2f\n",GNRL_ST.mfactor_th);
			if(GNRL_ST.calc_unique_flag)
				fprintf(fp,"                 Uniqueness >= %d\n\n",GNRL_ST.unique_th);
			else
				fprintf(fp,"                 Uniquness ignored\n\n");

			//dump top list
			if(motif_ids->size>MAX_FULL_MOTIF_LIST_SZ) {
				fprintf(fp,"\n\tTop Motifs List\n\n");
				fprintf(fp,"MOTIF\tCREAL\tCRAND\t\tCREAL\tCREAL\tUNIQ\tHITS\n");
				fprintf(fp,"ID\t[mili]\tSTATS\t\tZSCORE\tPVAL\t\n\n");
				for(l_id=list64_get_next(top_motif_ids,NULL); l_id!=NULL;
				l_id=list64_get_next(top_motif_ids,l_id)) {
					l_mtf=list64_get(res,l_id->val);
					mtf_res=(Motif_res*)l_mtf->p;

					if(mtf_res->conc_real_zscore!=UNDEFINED)
					fprintf(fp,"%lld\t%.3f\t%.3f+-%.3f\t%.2f\t%.3f\t%d\t%d\t\n",
						(uint64)mtf_res->id,
						mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
						mtf_res->conc_real_zscore,mtf_res->conc_real_pval,
						(int)mtf_res->unique_appear,(int)mtf_res->hits_num);
					else
						fprintf(fp,"%lld\t%.3f\t%.3f+-%.3f\t%.0f\t%.3f\t%d\t%d\t\n",
						(uint64)mtf_res->id,
						mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
						mtf_res->conc_real_zscore,mtf_res->conc_real_pval,
						(int)mtf_res->unique_appear,(int)mtf_res->hits_num);
					dump_motif_matrix(fp,mtf_res->id);
					fprintf(fp,"\n");
				}

				if(GNRL_ST.out_non_dangling_motifs) {
					fprintf(fp,"\n\n\tTotal No. of non-dangling motifs : %d\n\n",non_dangling_motif_ids->size);
					fprintf(fp,"MOTIF\tCREAL\tCRAND\t\tCREAL\tCREAL\tUNIQ\tHITS\n");
					fprintf(fp,"ID\t[mili]\tSTATS\t\tZSCORE\tPVAL\t\n\n");
					for(l_id=list64_get_next(non_dangling_motif_ids,NULL); l_id!=NULL;
					l_id=list64_get_next(non_dangling_motif_ids,l_id)) {
						l_mtf=list64_get(res,l_id->val);
						mtf_res=(Motif_res*)l_mtf->p;

						if(mtf_res->conc_real_zscore!=UNDEFINED)
						fprintf(fp,"%lld\t%.3f\t%.3f+-%.3f\t%.2f\t%.3f\t%d\t%d\t\n",
							(uint64)mtf_res->id,
							mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
							mtf_res->conc_real_zscore,mtf_res->conc_real_pval,
							(int)mtf_res->unique_appear,(int)mtf_res->hits_num);
						else
							fprintf(fp,"%lld\t%.3f\t%.3f+-%.3f\t%.0f\t%.3f\t%d\t%d\t\n",
							(uint64)mtf_res->id,
							mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
							mtf_res->conc_real_zscore,mtf_res->conc_real_pval,
							(int)mtf_res->unique_appear,(int)mtf_res->hits_num);
						dump_motif_matrix(fp,mtf_res->id);
						fprintf(fp,"\n");
					}
				} else {
					fprintf(fp,"\n\n\t(Total No. of non-dangling motifs : %d)\n\n",non_dangling_motif_ids->size);
				}
			}
			//dump all motifs
			fprintf(fp,"\n\n\tFull list includes %d motifs\n\n",motif_ids->size);
			fprintf(fp,"MOTIF\tCREAL\tCRAND\t\tCREAL\tCREAL\tUNIQ\tHITS\n");
			fprintf(fp,"ID\t[mili]\tSTATS\t\tZSCORE\tPVAL\t\n\n");
			for(l_id=list64_get_next(motif_ids,NULL); l_id!=NULL;
			l_id=list64_get_next(motif_ids,l_id)) {
				l_mtf=list64_get(res,l_id->val);
				mtf_res=(Motif_res*)l_mtf->p;

				if(mtf_res->conc_real_zscore!=UNDEFINED)
				fprintf(fp,"%lld\t%.3f\t%.3f+-%.3f\t%.2f\t%.3f\t%d\t%d\t\n",
					(uint64)mtf_res->id,
					mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
					mtf_res->conc_real_zscore,mtf_res->conc_real_pval,
					(int)mtf_res->unique_appear,(int)mtf_res->hits_num);
				else
					fprintf(fp,"%lld\t%.3f\t%.3f+-%.3f\t%.0f\t%.3f\t%d\t%d\t\n",
					(uint64)mtf_res->id,
					mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
					mtf_res->conc_real_zscore,mtf_res->conc_real_pval,
					(int)mtf_res->unique_appear,(int)mtf_res->hits_num);
				dump_motif_matrix(fp,mtf_res->id);
				fprintf(fp,"\n");
			}
		} else {
			//num of rand networks>=1000
			//consider pval
			fprintf(fp,"Criteria taken : Pval < %.3f\n",GNRL_ST.pval_th);
			fprintf(fp,"                 Mfactor > %.2f\n",GNRL_ST.mfactor_th);
			if(GNRL_ST.calc_unique_flag)
				fprintf(fp,"                 Uniqueness >= %d\n\n",GNRL_ST.unique_th);
			else
				fprintf(fp,"                 Uniquness ignored\n\n");

			fprintf(fp,"MOTIF\tCREAL\tCRAND\t\tCREAL\tCREAL\tUNIQ\n");
			fprintf(fp,"ID\t[mili]\tSTATS\t\tZSCORE\tPVAL\t\n");

			for(l_id=list64_get_next(motif_ids,NULL); l_id!=NULL;
			l_id=list64_get_next(motif_ids,l_id)) {
				l_mtf=list64_get(res,l_id->val);
				mtf_res=(Motif_res*)l_mtf->p;

				if(mtf_res->conc_real_zscore!=UNDEFINED)
				fprintf(fp,"%lld\t%.3f\t%.3f+-%.3f\t%.2f\t%.3f\t%d\t\n",
					(uint64)mtf_res->id,
					mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
					mtf_res->conc_real_zscore,mtf_res->conc_real_pval,
					(int)mtf_res->unique_appear);
				else
					fprintf(fp,"%lld\t%.3f\t%.3f+-%.3f\t%.0f\t%.3f\t%d\t\n",
					(uint64)mtf_res->id,
					mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
					mtf_res->conc_real_zscore,mtf_res->conc_real_pval,
					(int)mtf_res->unique_appear);
				dump_motif_matrix(fp,mtf_res->id);
				fprintf(fp,"\n");
			}
		}
		fprintf(fp,"\n");
	}




	//Dump full subgraphs list result to output file


	if (GNRL_ST.mtf_sz<=4) {
		//If motif size <=4 output all possible subgraphs list
		fprintf(fp,"\n\nFull list of subgraphs size %d ids:\n\n", GNRL_ST.mtf_sz);
		fprintf(fp,"\t( Total num of different subgraphs size %d is : %d )\n\n",GNRL_ST.mtf_sz,full_res->size);

		//deterministic - dump all id according to counts
		if(GNRL_ST.run_prob_app==FALSE) {
			fprintf(fp,"MOTIF\tNREAL\tNRAND\t\tNREAL\tNREAL\tCREAL\tUNIQ\n");
			fprintf(fp,"ID\t\tSTATS\t\tZSCORE\tPVAL\t[MILI]\t\n");
			for(l_mtf=list64_get_next(full_res,NULL); l_mtf!=NULL;l_mtf=list64_get_next(full_res,l_mtf)) {
				mtf_res=(Motif_res*)l_mtf->p;
				//eff = mtf_res->eff_arr[GNRL_ST.actual_n_eff];
				if(mtf_res->real_zscore!=UNDEFINED)
				fprintf(fp,"%lld\t%.0f\t%.1f+-%.1f\t%.2f\t%.3f\t%.2f\t%d\n",
					(uint64)mtf_res->id,
					mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
					mtf_res->real_zscore,mtf_res->real_pval,
					mtf_res->conc_real,(int)mtf_res->unique_appear);
				else
					fprintf(fp,"%lld\t%.0f\t%.1f+-%.1f\t%.0f\t%.3f\t%.2f\t%d\n",
					(uint64)mtf_res->id,
					mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
					mtf_res->real_zscore,mtf_res->real_pval,
						mtf_res->conc_real,(int)mtf_res->unique_appear);

				fprintf(fp,"\n");
			}
		//probabilistic
		}else{
			//print all id results according to concentrations
			fprintf(fp,"MOTIF\tCREAL\tCRAND\t\tCREAL\tCREAL\tUNIQ\tHITS\n");
			fprintf(fp,"ID\t[mili]\tSTATS\t\tZSCORE\tPVAL\t\t\n");

			for(l_mtf=list64_get_next(full_res,NULL); l_mtf!=NULL;
			l_mtf=list64_get_next(full_res,l_mtf)) {
				mtf_res=(Motif_res*)l_mtf->p;
				if(mtf_res->conc_real_zscore!=UNDEFINED)
				fprintf(fp,"%lld\t%.3f\t%.3f+-%.3f\t%.2f\t%.3f\t%d\t",
					(uint64)mtf_res->id,
					mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
					mtf_res->conc_real_zscore,mtf_res->conc_real_pval,(int)mtf_res->unique_appear);
				else
					fprintf(fp,"%lld\t%.3f\t%.3f+-%.3f\t%.0f\t%.3f\t%d\t",
					(uint64)mtf_res->id,
					mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
						mtf_res->conc_real_zscore,mtf_res->conc_real_pval,(int)mtf_res->unique_appear);

				fprintf(fp,"%d\n",mtf_res->hits_num);
				fprintf(fp,"\n");
			}
		}
	} else {
		//if motif size >=5 output only full list of subgraph exists in the real network
		fprintf(fp,"\n\nFull list of subgraphs size %d found in the network:\n\n", GNRL_ST.mtf_sz);
		fprintf(fp,"\t( Total num of different subgraphs size %d found in the real network : %d )\n\n",GNRL_ST.mtf_sz,res->size);

		//deterministic - dump all id according to counts
		if(GNRL_ST.run_prob_app==FALSE) {
			fprintf(fp,"MOTIF\tNREAL\tNRAND\t\tNREAL\tNREAL\tCREAL\tUNIQ\n");
			fprintf(fp,"ID\t\tSTATS\t\tZSCORE\tPVAL\t[MILI]\t\n");
			for(l_mtf=list64_get_next(res,NULL); l_mtf!=NULL;l_mtf=list64_get_next(res,l_mtf)) {
				mtf_res=(Motif_res*)l_mtf->p;
				if(mtf_res->real_zscore!=UNDEFINED)
				fprintf(fp,"%lld\t%.0f\t%.1f+-%.1f\t%.2f\t%.3f\t%.2f\t%d\n",
					(uint64)mtf_res->id,
					mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
					mtf_res->real_zscore,mtf_res->real_pval,
					mtf_res->conc_real,(int)mtf_res->unique_appear);
				else
					fprintf(fp,"%lld\t%.0f\t%.1f+-%.1f\t%.0f\t%.3f\t%.2f\t%d\n",
					(uint64)mtf_res->id,
					mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
					mtf_res->real_zscore,mtf_res->real_pval,
						mtf_res->conc_real,(int)mtf_res->unique_appear);

				fprintf(fp,"\n");
			}
		//probabilistic
		}else{
			//print all id results according to concentrations
			fprintf(fp,"MOTIF\tCREAL\tCRAND\t\tCREAL\tCREAL\tUNIQ\tHITS\n");
			fprintf(fp,"ID\t[mili]\tSTATS\t\tZSCORE\tPVAL\t\t\n");

			for(l_mtf=list64_get_next(res,NULL); l_mtf!=NULL;
			l_mtf=list64_get_next(res,l_mtf)) {
				mtf_res=(Motif_res*)l_mtf->p;

				if(mtf_res->conc_real_zscore!=UNDEFINED)
				fprintf(fp,"%lld\t%.3f\t%.3f+-%.3f\t%.2f\t%.3f\t%d\t",
					(uint64)mtf_res->id,
					mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
					mtf_res->conc_real_zscore,mtf_res->conc_real_pval,(int)mtf_res->unique_appear);
				else
					fprintf(fp,"%lld\t%.3f\t%.3f+-%.3f\t%.0f\t%.3f\t%d\t",
					(uint64)mtf_res->id,
					mtf_res->conc_real,mtf_res->conc_rand_mean, mtf_res->conc_rand_std_dev,
						mtf_res->conc_real_zscore,mtf_res->conc_real_pval,(int)mtf_res->unique_appear);
				fprintf(fp,"%d\n",mtf_res->hits_num);
				fprintf(fp,"\n");
			}
		}
	}
	fprintf(fp,"\n");
	if( motif_ids->size>=MAX_FULL_MOTIF_LIST_SZ) {
		list64_free_mem(top_motif_ids);
	}
	list64_free_mem(non_dangling_motif_ids);
	list64_free_mem(motif_ids);
	fflush(fp);
}


//output all ids data matrix short format- good for matlab use
//in counts units
void
dump_all_ids_N_data_matrix(FILE *fp, list64* full_res)
{
	list64_item *l_mtf;
	Motif_res *mtf_res;

	for(l_mtf=list64_get_next(full_res,NULL); l_mtf!=NULL;l_mtf=list64_get_next(full_res,l_mtf)) {
		mtf_res=(Motif_res*)l_mtf->p;
		fprintf(fp,"%lld %.0f %.4f %.4f %.4f %.4f %d %.4f\n",
			(uint64)mtf_res->id,
			mtf_res->real_count,mtf_res->rand_mean, mtf_res->rand_std_dev,
			mtf_res->real_zscore,mtf_res->real_pval,(int)mtf_res->unique_appear,
			(double)mtf_res->conc_real);
	}
}


//output all ids data matrix short format- good for matlab use
//in conc units
void
dump_all_ids_C_data_matrix(FILE *fp,list64* full_res)
{
	list64_item *l_mtf;
	Motif_res *mtf_res;

	for(l_mtf=list64_get_next(full_res,NULL); l_mtf!=NULL;l_mtf=list64_get_next(full_res,l_mtf)) {
		mtf_res=(Motif_res*)l_mtf->p;
		fprintf(fp,"%lld\t%.5f\t%.5f\t%.5f\t%.2f\t%.3f\t%d\n",
			(uint64)mtf_res->id,
			mtf_res->conc_real,
			mtf_res->conc_rand_mean,
			mtf_res->conc_rand_std_dev,
			mtf_res->conc_real_zscore,
			mtf_res->conc_real_pval,
			(int)mtf_res->unique_appear);
	}

}




void
dump_rand_all_to_mat_file(FILE *fp,int rnd_net_num,list64 *full_res_all)
{
	int i;
	list64_item *l_i;
	//list *net_res;
	Motif_res *mtf_res;

	// dump all rand results
	for( l_i=list64_get_next(full_res_all,NULL);l_i!=NULL;l_i=list64_get_next(full_res_all,l_i)) {
		mtf_res=(Motif_res*)l_i->p;
		for(i=0; i<rnd_net_num; i++) {
			fprintf(fp,"%d\t", mtf_res->all_rand_counts[i]);
		}
		fprintf(fp,"\n");
	}

}


// dump run time measurement results
//arguments:
//	fp - file pointer to dump to
// text to print
// runtime structure to dump time measuremnt from
void
dump_time_measure(FILE* fp,char *text, Runtime *runtime)
{
	double run_time = runtime->elapsed;

	if (runtime->elapsed > 3600) {
		//display message in hours
		run_time /=(double)3600;
		fprintf(fp,"\n (%s %6.2f hours.)\n", text, run_time);
	}
	else if (run_time > 60) {
		//display message in minutes
		run_time /=(double)60;
		fprintf(fp,"\n (%s %6.2f minutes).\n", text, run_time);
	}
	else {
		//display message in seconds
		fprintf(fp,"\n (%s %6.1f seconds.)\n", text, run_time );
	}
}


void
print_public_help()
{

	printf("Usage : mfinder <Network input file name> -s <motif size> -r <no. of randon networks> [-f <output file name>] [more flags]");
	printf("\n\n");
	printf("\t-s <motif size>  :Motif size to search\n");
	printf("\t-r <rand net num> :Number of random networks to generate\n");
	printf("\t-f <output file name>  : Output file name\n");
	printf("\t-nd : Input network is a non-directed network.\n");
	printf("\t-p <num of samples>: run with Sampling method,\n");
	printf("\t-omem : output members list of all subgraphs\n");
	printf("\t-h : help\n");

	printf("\n\tAdditional flags:\n");

	printf("\n\tMotif criteria flags:\n");
	printf("\t-m <value> : mfactor threshold to use when calculating motifs\n");
	printf("\t-z <value> : Z-score threshold to use when calculating motifs\n");
	printf("\t-u : Uniqueness threshold\n");
	printf("\t-nu : Dont count uniqueness and ignore uniqueness threshold\n");

	printf("\n\tRandom networks randomization flags:\n\n");
	printf("\t-rs : use stubs method for generating random networks\n");
	printf("\t-rclust : Preserve clustering sequence in random networks\n");
	printf("\t-met :Use Metropolis algorithm to conserve triad-census\n\t\tin random networks\n");
	printf("\t\t(for s>3; Default : Do not use Metropolis)\n");
	printf("\t-t0 <(default 0.001)> :Initial temperature (-met option)\n");
	printf("\t-iter <(default 2)> :controls how many steps to perform (-met option) \n");
	printf("\t-eth  <(default 0.005)> : energy threshhold (-met option)\n");
	printf("\t-rgrass <colony size>: generate random networks using grassberger \n\t\talgorithm\n");
	printf("\t-rgrass_max_sz <max ratio>: Limit maximal colony size ratio\n");
	printf("\t-rdm: don't conserve mutual edges in random networks\n");
	printf("\t-rcl <layers num><size1 size2 ..sizem>: conserve layers in random\n\t\t networks\n");
	printf("\t-nsr : Global Switches number when generating Random networks.\n");
	printf("\n");

	printf("\n\tOutput files flags:\n\n");
	printf("\t-oi : output intermediate output file. Defualt :No: \n");
	printf("\t-ospmem <subgraph id>: output members list of a specific subgraphs only\n");
	printf("\t-maxmem <list length>: limit length of members list to 'list length'.\n\t\tDefualt: 1000\n");
	printf("\t-omat : output matrix format file ('__MAT.txt')\n");
	printf("\t-omet : output metropolis log\n");
	printf("\t-olog : output general log file\n");
	printf("\t-orall : output matrix format of appearances in each random network\n");
	printf("\t-ornet : output random networks files\n");
	printf("\t-otop <no. of top motifs> : No. of top motifs to show\n");
	printf("\t-onodangl: output a list of all nan-dangling detected motifs\n");

	printf("\n\tOther flags\n\n");
	printf("\t-ts : <target,source,weight> Old format of input network file\n");
	printf("\t-q :Quiet mode - No output to the screen\n");
	printf("\t-dd : Don't die mode. Wait to user action before terminating the\n\t\tprogram\n");
	printf("\t-pold <num of samples>: run sampling method old version\n");
	printf("\t-nor : Dont search Real network. Defualt :No: \n");
	printf("\t-cr : calculate roles statistics\n");
}

void
print_private_help()
{
	printf("Usage : mfinder <Network input file name> -s <motif size> -r <rand net num> [-f <output file name>] [-u]");
	printf("\n");
	printf("\t-s <motif size>  :motif size to search\n");
	printf("\t-r <rand net num> :number of random network to generate\n");
	printf("\t-f <output file name>  : Output file name\n");
	printf("\t-ts : <target, source> Old format of input network file\n");
	printf("\t-se : Input network includes self edges.\n");
	printf("\t-w :  Input network includes weights.\n");
	printf("\n");
	printf("\t-u : uniqueness threshold\n");
	printf("\t-nu : Dont count uniqueness and ignore uniqueness threshold\n");
	printf("\t-m <value> : mfactor threshold to use when calculating motifs\n");
	printf("\t-z <value> : Z-score threshold to use when calculating motifs\n");
	printf("\t-cr : calculate roles statistics\n");
	printf("\t-rs : use stubs method for generating random networks\n");
	printf("\t-rclust : Preserve clustering sequence in random networks\n");
	printf("\n");
	printf("\t-pold <num of samples>: run probabilistic approach old version\n");
	printf("\t-p <num of samples>: run probabilistic approach, efficient version\n");
	printf("\t-met :Use Metropolis mode (Default : Do not use Metropolis)\n");
	printf("\t-t0 <(default 0.001)> :Initial temperature (-met option)\n");
	printf("\t-iter <(default 2)> :controls how many steps to perform (-met option) \n");
	printf("\t-eth  <(default 0.005)> : energy threshhold (-met option)\n");
	printf("\n");
	printf("\t-nor : Dont search Real network. Defualt :No: \n");
	printf("\t-nsr : Global Switches number when generating Random networks (Consider unsuccessfull switches).\n");
	printf("\n");
	printf("\t-q :Quiet mode - No output to the screen\n");
	printf("\t-dd : Don't die mode. wait to user action before terminating the program\n");
	printf("\t-oi : output intermediate output file. Defualt :No: \n");
	printf("\t-ol : output matlab matrix long format\n");
	printf("\t-oc : output matlab short format with Concentrations and with Counts\n");
	printf("\t-omem : output members list of all subgraphs\n");
	printf("\t-ospmem <subgraph id>: output members list of a specific subgraphs only\n");
	printf("\t-maxmem <list length>: limit length of members list to 'list length'. Defualt: 1000\n");
	printf("\t-omat : output matlab short format file\n");
	printf("\t-omet : output metropolis log\n");
	printf("\t-olog : output general log file\n");
	printf("\t-orall : output matlab format all random results\n");
	printf("\t-ornet : output random networks files\n");
	printf("\t-oclust : output clustering sequences of real and random nets\n");
	printf("\t-otop <no. of top motifs> : No. of top motifs to show\n");
	printf("\t-onodangl: output a list of all nan-dangling detected motifs\n");
	printf("\t-rgrass <colony size>: generate random network using grassberger algorithm\n");
	printf("\t-rgrass_max_sz <max ratio>: Limit maximal colony size ratio\n");
	printf("\t-rdm: don't conserve mutual edges in random networks\n");
	printf("\t-rcl <layers num><size1 size2 ..sizem>: conserve layers in random networks\n");
	printf("\t-h : help\n");
	printf("\t-hh : hidden help\n");
}



int
output_results(list64 *final_res, list64*final_res_all)
{
	int rc=RC_OK;

	dump_final_res(GNRL_ST.out_fp,final_res,final_res_all,GNRL_ST.rnd_net_num);

	//if not finished generating  random networks , retunr here
	//finished all random networks generating
	if(GNRL_ST.out_s_mat_flag==TRUE) {
		if(GNRL_ST.out_s_c_mat_flag==TRUE){
			if(GNRL_ST.mtf_sz<=4)
				dump_all_ids_C_data_matrix(GNRL_ST.mat_s_fp,final_res_all);
			else
				//then final_res_all is empty , so send with final_res instead
				dump_all_ids_C_data_matrix(GNRL_ST.mat_s_fp,final_res);
		}else{
			if(GNRL_ST.mtf_sz<=4)
				dump_all_ids_N_data_matrix(GNRL_ST.mat_s_fp,final_res_all);
			else
				dump_all_ids_N_data_matrix(GNRL_ST.mat_s_fp,final_res);
		}
	}
	if(GNRL_ST.out_rand_mat_flag==TRUE)
			dump_rand_all_to_mat_file(GNRL_ST.mat_rand_fp,GNRL_ST.rnd_net_num,final_res_all);
	if(GNRL_ST.out_members==TRUE)
		dump_members_lists(GNRL_ST.members_fp,final_res);
	if(GNRL_ST.out_roles==TRUE)
		dump_roles_results(GNRL_ST.roles_fp,final_res_all,GNRL_ST.rnd_net_num);
	if(GNRL_ST.out_clustering==TRUE)
		dump_clustering_series(G_N,GNRL_ST.clust_fp);
	//output time measure
	//random network elapsed time should be diveided by random networks num
	if(GNRL_ST.rnd_net_num>0)
		GNRL_ST.rand_net_time.elapsed /= (double)GNRL_ST.rnd_net_num;

	if(GNRL_ST.quiet_mode==FALSE) {
		dump_time_measure(stdout, "Application total runtime was:", &GNRL_ST.total_time);
		dump_time_measure(stdout, "Real network processing runtime was:", &GNRL_ST.real_net_time);
		dump_time_measure(stdout, "Single Random network processing runtime was:", &GNRL_ST.rand_net_time);
	}

	dump_time_measure(GNRL_ST.out_fp, "Application total runtime was:", &GNRL_ST.total_time);
	dump_time_measure(GNRL_ST.out_fp, "Real network processing runtime was:", &GNRL_ST.real_net_time);
	dump_time_measure(GNRL_ST.out_fp, "Single Random network processing runtime was:", &GNRL_ST.rand_net_time);

	fclose(GNRL_ST.out_fp);
	if(GNRL_ST.out_log_flag==TRUE)
		fclose(GNRL_ST.log_fp);
	if(GNRL_ST.out_s_mat_flag==TRUE)
		fclose(GNRL_ST.mat_s_fp);
	if(GNRL_ST.out_metrop_mat_flag==TRUE)
		fclose(GNRL_ST.mat_metrop_fp);
	if(GNRL_ST.out_l_mat_flag==TRUE)
		fclose(GNRL_ST.mat_l_fp);
	if(GNRL_ST.out_members==TRUE)
		fclose(GNRL_ST.members_fp);
	if(GNRL_ST.out_roles==TRUE)
		fclose(GNRL_ST.roles_fp);
	if(GNRL_ST.out_clustering==TRUE)
		fclose(GNRL_ST.clust_fp);
	if(GNRL_ST.out_rand_mat_flag==TRUE)
		fclose(GNRL_ST.mat_rand_fp);

	if(GNRL_ST.quiet_mode==FALSE) {
		if(GNRL_ST.out_s_mat_flag==TRUE)
			printf("\nMatrix output file %s was generated\n",GNRL_ST.mat_s_fname);
		if(GNRL_ST.out_l_mat_flag==TRUE)
			printf("Long Matrix output file %s was generated\n",GNRL_ST.mat_l_fname);
		printf("\nOutput File %s was generated\n",GNRL_ST.out_fname);
	}
	return rc;
}


void
output_network_to_text_file(Network *N, char *fname)
{
	FILE *fp;
	int i,j;
	int val;

	fp=fopen(fname,"wt");
	if(fp==NULL){
		printf("Error: Check output file path\n");
		printf("%s\n",fname);
		//return(0);
		at_exit(-1);
	}
	for (i=1;i<=N->vertices_num;i++){
		for (j=1;j<=N->vertices_num;j++){
			if( (val=MatGet(N->mat,i,j))!=0) {
				if(GNRL_ST.input_net_format==SRC_TRG_FORMAT)
					fprintf(fp,"%d %d %d\n",i,j,val);
				else
					fprintf(fp,"%d %d %d\n",j,i,val);
			}
		}
	}


	fclose(fp);
}





