/************************************************************************
*
*  File name: results.c
*
*  Description: functions realted to results calculations
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "common.h"
#include "list.h"
#include "mat.h"
#include "hash.h"
#include "motif_ids.h"
#include "metropolis.h"
#include "role.h"
#include "random.h"





/*************************** Global variables ****************************/


/******************************* Externs *********************************/
extern int DEBUG_LEVEL;

extern char *input_network_fname;
extern Network *G_N;
extern Gnrl_st GNRL_ST;
//result table
extern Res_tbl RES_TBL;
//Gneralized motifs res table
extern Res_tbl GMTF_RES_TBL,GMTF_CMPLX_RES_TBL;

extern FILE *out_fp,*log_fp;
extern time_t start_time, end_time;

extern Hash *Role_hash;

extern list64 *final_res, *final_res_all;
extern list64 *res_sub_motif;

extern int *Perms;



/******************************* Functions *******************************/




void
init_res_tbl(Res_tbl *res)
{
	int i;
	list64_init(&res->real);
	res->rand_arr=(list64**)calloc(GNRL_ST.rnd_net_num,sizeof(list64*));
	for(i=0; i<GNRL_ST.rnd_net_num; i++)
		list64_init((list64**)&res->rand_arr[i]);
}

void
res_tbl_mem_free_single(list64 *L)
{
	list64_item *item;
	Motif*	mtf;

	for(item=list64_get_next(L,NULL);item!=NULL;item=list64_get_next(L,item)) {
		if (item->p != NULL)	{
			mtf = (Motif*)item->p;
			if (mtf->members != NULL)
				list_free_mem(mtf->members);
			if (mtf->all_members != NULL)
				list_free_mem(mtf->all_members);
		}
	}
	list64_free_mem(L);
}

void
res_tbl_mem_free(Res_tbl *res)
{
	int i;

// handling results of random networks
	for(i=0;i<GNRL_ST.rnd_net_num;i++)
		res_tbl_mem_free_single (res->rand_arr[i]);
	free(res->rand_arr);

// handling results of the real network
	res_tbl_mem_free_single (res->real);
}

void
final_res_free(list64 *motif_res_lis)
{
	if(motif_res_lis != NULL) {
		list64_free_mem(motif_res_lis);
	}
}


void
calc_stat(double*arr, int sz, double* mean_p, double *stddev_p)
{
	int i;
	double sum=0;
	double mean=0;
	double stddev=0;
	double sigma=0;
	double power=0;

	//calc mean val
	for(i=0; i<sz; i++)
		sum += arr[i];
	if(sz > 0)
		mean = (double)(sum / (double)sz);
	else
		mean=0;
	//calc stddev
	if(sz >1) {
		for(i=0; i<sz; i++) {
			if( (arr[i]-mean) == 0)
				power = 0;
			else
				power=pow((double)((double)arr[i] - (double)mean), (double)2);
			sigma += power/(double)(sz);
		}
	}
	stddev=(double)sqrt(sigma);

	*mean_p=mean;
	*stddev_p=stddev;
}



void
join_subgraphs_res(list64 **res_p, int mtf_sz,int rand_net_indx)
{
	list_item *l_members;
	list64_item *l_res_id, *l_res,*l_iso,*l_iso_rep;
	int64 rep_id;
	int64 curr_id;
	list64 *iso_id_list;
	list *members_list;
	list64 *final_res,*res=*res_p;
	Motif *rep_mtf;
	int i,unique;
	int new_members;
	int *mtf_vrtx_arr, *new_vrtx_arr;
	Member *mtf_members_arr, *new_members_arr,*new_members_arr_;
	double total_prob_counts=0;
	double total_motifs_counts=0;
	int curr_perm_index,rep_id_perm_index;
	int *curr_perm,*rep_id_perm;

	list64_init(&final_res);

	//pass through res list and for each id look for the isomorphism ids
	while((l_res_id=list64_get_next(res,NULL)) != NULL) {
			iso_id_list=calc_mtf_id_iso((int64)l_res_id->val, mtf_sz);
			//find all other iso id, remove them from res list
			//and sum all iso id counts
			//insert new item with iso_id (first id in iso_id_list) to final_res
			l_iso_rep=list64_get_next(iso_id_list,NULL);
			rep_id=(int64)l_iso_rep->val;
			rep_mtf=(Motif*)calloc(1,sizeof(Motif));
			
			rep_mtf->id = rep_id;
			list_init(&rep_mtf->members);
			list_init(&rep_mtf->all_members);
			//info of the permuations of the rep id
			//needed for mapping the members
			rep_id_perm_index=*(int*)l_iso_rep->p;
			rep_id_perm=get_perm(mtf_sz,rep_id_perm_index);

			//pass through iso_id_list
			for(l_iso=list64_get_next(iso_id_list,NULL); l_iso!=NULL;
				l_iso=list64_get_next(iso_id_list,l_iso)) {
					curr_id=(int64)l_iso->val;
					//get the permutaqiton in order to sort the members in their corresponding order (defined by the ID)
					curr_perm_index=*(int*)l_iso->p;
					//for(i=0;i<mtf_sz;i++)
						//perm[i]=Perms[mtf_sz][curr_perm_index][i];
					curr_perm=get_perm(mtf_sz,curr_perm_index);
					//find element in res
					if( (l_res=list64_get(res,(int64)curr_id)) != NULL) {
						//add to score of rep id
						rep_mtf->count+=(double)((Motif*)(l_res->p))->count;
						rep_mtf->hits+=((Motif*)(l_res->p))->hits;
						rep_mtf->prob_count+=((Motif*)(l_res->p))->prob_count;
						//if calc unique
						if(GNRL_ST.calc_unique_flag == TRUE) {
							//merge list of unique memebers of curr_id to rep_id
							if( (members_list=((Motif*)(l_res->p))->members) != NULL) {

								for(l_members = list_get_next(members_list,NULL); l_members!=NULL;
								l_members = list_get_next(members_list,l_members)) {
									//pass through list , if this vertices group is unique (no common vertices with
									//same motif id that already found) then insert to list
									mtf_vrtx_arr=(int*)l_members->p;
									unique = check_if_unique(rep_mtf,mtf_vrtx_arr,mtf_sz);
									if(unique==TRUE) {
										new_vrtx_arr=(int*)calloc(mtf_sz,sizeof(int));
										
										//copy it cause later i free all members of l_res motif
										for(i=0;i<mtf_sz;i++)
											new_vrtx_arr[i]=mtf_vrtx_arr[i];
										list_insert(rep_mtf->members,rep_mtf->members->size+1,(void*)new_vrtx_arr);
									}
								}
								//free members list from curr_id
								list_free_members_list(members_list);
							}
						}
						if(GNRL_ST.list_members == TRUE) {
							//merge list of members of curr_id to rep_id
							if( (members_list=((Motif*)(l_res->p))->all_members) != NULL) {

								for(l_members = list_get_next(members_list,NULL); l_members!=NULL;
								l_members = list_get_next(members_list,l_members)) {
									//pass through list , if this vertices group is unique (no common vertices with
									//same motif id that already found) then insert to list
									mtf_members_arr=(Member*)l_members->p;
									new_members=check_if_new_members(rep_mtf,mtf_members_arr,mtf_sz);
									if(new_members==TRUE) {
										new_members_arr_=(Member*)calloc(mtf_sz,sizeof(Member));
										new_members_arr=(Member*)calloc(mtf_sz,sizeof(Member));
										//copy it cause later i free all members of l_res motif

										//copy in the right order to fit the motif ID
										//for that I use the permutaiton of the specific id
										//this is done in tw steps:
										//first map to zero perm
										for(i=0;i<mtf_sz;i++){

											new_members_arr_[curr_perm[i]].node=mtf_members_arr[i].node;
											new_members_arr_[curr_perm[i]].role=mtf_members_arr[i].role;

										}
										for(i=0;i<mtf_sz;i++){
										//second map by rep_id_perm
											new_members_arr[i].node=new_members_arr_[rep_id_perm[i]].node;
											new_members_arr[i].role=new_members_arr_[rep_id_perm[i]].role;
										}
										free(new_members_arr_);
										list_insert(rep_mtf->all_members,1,(void*)new_members_arr);
									}
								}
								//free members list from curr_id
								list_free_all_members_list(members_list);
							}
						}
						members_list=((Motif*)(l_res->p))->members;
						list_free_mem(members_list);
						members_list=((Motif*)(l_res->p))->all_members;
						list_free_mem(members_list);
						//remove curr_id item from res list
						list64_delete(res,(int64)curr_id);
					}
			}
			list64_free_mem(iso_id_list);
			list64_insert(final_res,(int64)rep_id,(void*)rep_mtf);
	}

	//if probabilistic alg then here is the place to renormalize
	//the prob_count field
	if(GNRL_ST.run_prob_app==TRUE) {
		//calc total_prob_count
		for(l_res_id=list64_get_next(final_res,NULL); l_res_id!=NULL;
			l_res_id=list64_get_next(final_res,l_res_id)) {
				rep_mtf=(Motif*)l_res_id->p;
				total_prob_counts+=rep_mtf->prob_count;
		}
		//calc count for each motif
		//calc = prob_count/toal_prob_count
		for(l_res_id=list64_get_next(final_res,NULL); l_res_id!=NULL;
			l_res_id=list64_get_next(final_res,l_res_id)) {
				rep_mtf=(Motif*)l_res_id->p;
				rep_mtf->count=(double)(((double)rep_mtf->prob_count/(double)total_prob_counts)*(double)GNRL_ST.prob_total_samples_num[rand_net_indx]);
				rep_mtf->conc=1000*(double)rep_mtf->count/(double)GNRL_ST.prob_total_samples_num[rand_net_indx];
		}
	} else {
		//deterministic alg - here is the place to caclulate the concentrations

		//calc total motifs count
		for(l_res_id=list64_get_next(final_res,NULL); l_res_id!=NULL;
			l_res_id=list64_get_next(final_res,l_res_id)) {
				rep_mtf=(Motif*)l_res_id->p;
				total_motifs_counts+=rep_mtf->count;
		}
		//calc concentrations
		for(l_res_id=list64_get_next(final_res,NULL); l_res_id!=NULL;
			l_res_id=list64_get_next(final_res,l_res_id)) {
				rep_mtf=(Motif*)l_res_id->p;
				rep_mtf->conc=1000*(double)rep_mtf->count/(double)total_motifs_counts;
		}
	}

	//free res list and make res point to final res list
	res_tbl_mem_free_single(res);
	//list64_free_mem(res);
	*res_p=final_res;

}


void
update_global_res_tbl(list64 *global_res,list64 *curr_res,int rand_net_indx)
{
	list64_item *glb_l,*curr_l;
	Motif *glb_mtf,*mtf;
	list *members_list;
	list_item *l_members;
	int *mtf_vrtx_arr,*new_vrtx_arr;
	int unique;
	Member *mtf_members_arr,*new_members_arr;
	int new_members,i;
	double total_motifs_counts=0;
	double total_prob_counts=0;

	for(curr_l=list64_get_next(curr_res,NULL);curr_l!=NULL;curr_l=list64_get_next(curr_res,curr_l)){
		mtf=(Motif*)curr_l->p;
		if( (glb_l=list64_get(global_res,mtf->id))==NULL) {
			//New subgraph
			glb_mtf=(Motif*)calloc(1,sizeof(Motif));
			glb_mtf->id = mtf->id;
			list_init(&glb_mtf->members);
			list_init(&glb_mtf->all_members);

			glb_mtf->count=mtf->count;
			glb_mtf->hits=mtf->hits;
			glb_mtf->prob_count=mtf->prob_count;
			glb_mtf->numberOfSelfEdges=mtf->numberOfSelfEdges;
			glb_mtf->conv_grade=mtf->conv_grade;
			//copy members list and all members list
			//if calc unique
			if(GNRL_ST.calc_unique_flag == TRUE) {
				//merge list of unique memebers of curr_id to rep_id
				if( (members_list=mtf->members) != NULL) {
					for(l_members = list_get_next(members_list,NULL); l_members!=NULL;
					l_members = list_get_next(members_list,l_members)) {
						//pass through list , if this vertices group is unique (no common vertices with
						//same motif id that already found) then insert to list
						mtf_vrtx_arr=(int*)l_members->p;
						new_vrtx_arr=(int*)calloc(GNRL_ST.mtf_sz,sizeof(int));
						
						//copy it cause later i free all members of l_res motif
						for(i=0;i<GNRL_ST.mtf_sz;i++)
							new_vrtx_arr[i]=mtf_vrtx_arr[i];
						list_insert(glb_mtf->members,mtf->members->size+1,(void*)new_vrtx_arr);
					}
					//free members list from curr_id
					list_free_members_list(members_list);
					list_free_mem(members_list);
				}
			}
			if(GNRL_ST.list_members == TRUE) {
				//merge list of members of curr_id to rep_id
				if( (members_list=mtf->all_members) != NULL) {
					for(l_members = list_get_next(members_list,NULL); l_members!=NULL;
					l_members = list_get_next(members_list,l_members)) {
						//pass through list , if this vertices group is unique (no common vertices with
						//same motif id that already found) then insert to list
						mtf_members_arr=(Member*)l_members->p;
						new_members_arr=(Member*)calloc(GNRL_ST.mtf_sz,sizeof(Member));

						//copy it cause later i free all members of l_res motif
						for(i=0;i<GNRL_ST.mtf_sz;i++){
							new_members_arr[i].node=mtf_members_arr[i].node;
							new_members_arr[i].role=mtf_members_arr[i].role;
						}
						list_insert(glb_mtf->all_members,1,(void*)new_members_arr);
					}
					//free members list from curr_id
					list_free_all_members_list(members_list);
					list_free_mem(members_list);
				}
			}
			list64_insert(global_res,glb_mtf->id,(void*)glb_mtf);
		}else{
			//Subgrpah already exist
			glb_mtf=(Motif*)glb_l->p;
			glb_mtf->count+=mtf->count;
			glb_mtf->hits+=mtf->hits;
			glb_mtf->prob_count+=mtf->prob_count;
			glb_mtf->numberOfSelfEdges=mtf->numberOfSelfEdges;
			glb_mtf->conv_grade=mtf->conv_grade;

			//if calc unique
			if(GNRL_ST.calc_unique_flag == TRUE) {
				//merge list of unique memebers of curr_id to rep_id
				if( (members_list=mtf->members) != NULL) {
					for(l_members = list_get_next(members_list,NULL); l_members!=NULL;
					l_members = list_get_next(members_list,l_members)) {
						//pass through list , if this vertices group is unique (no common vertices with
						//same motif id that already found) then insert to list
						mtf_vrtx_arr=(int*)l_members->p;
						unique = check_if_unique(glb_mtf,mtf_vrtx_arr,GNRL_ST.mtf_sz);
						if(unique==TRUE) {
							new_vrtx_arr=(int*)calloc(GNRL_ST.mtf_sz,sizeof(int));
							
							//copy it cause later i free all members of l_res motif
							for(i=0;i<GNRL_ST.mtf_sz;i++)
								new_vrtx_arr[i]=mtf_vrtx_arr[i];
							list_insert(glb_mtf->members,mtf->members->size+1,(void*)new_vrtx_arr);
						}
					}
					//free members list from curr_id
					list_free_members_list(members_list);
					list_free_mem(members_list);
				}
			}
			if(GNRL_ST.list_members == TRUE) {
				//merge list of members of curr_id to rep_id
				if( (members_list=mtf->all_members) != NULL) {
					for(l_members = list_get_next(members_list,NULL); l_members!=NULL;
					l_members = list_get_next(members_list,l_members)) {
						//pass through list , if this vertices group is unique (no common vertices with
						//same motif id that already found) then insert to list
						mtf_members_arr=(Member*)l_members->p;
						new_members=check_if_new_members(glb_mtf,mtf_members_arr,GNRL_ST.mtf_sz);
						if(new_members==TRUE) {
							new_members_arr=(Member*)calloc(GNRL_ST.mtf_sz,sizeof(Member));

							//copy it cause later i free all members of l_res motif
							for(i=0;i<GNRL_ST.mtf_sz;i++){
								new_members_arr[i].node=mtf_members_arr[i].node;
								new_members_arr[i].role=mtf_members_arr[i].role;
							}
							list_insert(glb_mtf->all_members,1,(void*)new_members_arr);
						}
					}
					//free members list from curr_id
					list_free_all_members_list(members_list);
					list_free_mem(members_list);
				}
			}
		}
	}
	//if probabilistic alg then here is the place to renormalize
	//the prob_count field
	if(GNRL_ST.run_prob_app==TRUE) {
		//calc total_prob_count
		for(glb_l=list64_get_next(global_res,NULL); glb_l!=NULL;
		glb_l=list64_get_next(global_res,glb_l)) {
				glb_mtf=(Motif*)glb_l->p;
				total_prob_counts+=glb_mtf->prob_count;
		}
		//calc count for each motif
		for(glb_l=list64_get_next(global_res,NULL); glb_l!=NULL;
		glb_l=list64_get_next(global_res,glb_l)) {
			//calc = prob_count/toal_prob_count
			glb_mtf=(Motif*)glb_l->p;
			glb_mtf->count=(double)(((double)glb_mtf->prob_count/(double)total_prob_counts)*(double)GNRL_ST.prob_total_samples_num[rand_net_indx]);
			glb_mtf->conc=1000*(double)glb_mtf->count/(double)GNRL_ST.prob_total_samples_num[rand_net_indx];
		}
	} else {
		//deterministic alg - here is the place to caclulate the concentrations
		for(glb_l=list64_get_next(global_res,NULL); glb_l!=NULL;
		glb_l=list64_get_next(global_res,glb_l)) {
			glb_mtf=(Motif*)glb_l->p;
			total_motifs_counts+=glb_mtf->count;
		}
		//calc concentrations
		for(glb_l=list64_get_next(global_res,NULL); glb_l!=NULL;
		glb_l=list64_get_next(global_res,glb_l)) {
			//calc = prob_count/toal_prob_count
			glb_mtf=(Motif*)glb_l->p;
			glb_mtf->conc=1000*(double)glb_mtf->count/(double)total_motifs_counts;
		}
	}
}



/****************** statistics calculation for results ******************/


void
calc_deviation(Motif_res *mtf_res, Res_tbl *res_tbl, int rnd_net_num)
{
	int i;
	list64_item *l_mtf;
	Motif *mtf;
	double counts_sum=0,x;
	double conc_sum=0,y;
	double sigma=0;
	double power;

	// calc counts mean val
	for(i=0; i<rnd_net_num; i++) {
		if( (l_mtf=list64_get(res_tbl->rand_arr[i],(int64)mtf_res->id)) !=NULL) {
			mtf=(Motif*)l_mtf->p;
			counts_sum += mtf->count;
			if(GNRL_ST.out_rand_mat_flag)
				mtf_res->all_rand_counts[i]=(int)mtf->count;
		}
	}
	if(rnd_net_num > 0)
		mtf_res->rand_mean = (double)counts_sum / (double)(rnd_net_num);
	else
		mtf_res->rand_mean=0;
	// calc concentration mean val
	for(i=0; i<rnd_net_num; i++) {
		if( (l_mtf=list64_get(res_tbl->rand_arr[i],(int64)mtf_res->id)) !=NULL) {
			mtf=(Motif*)l_mtf->p;
			conc_sum += mtf->conc;
		}
	}
	if(rnd_net_num > 0)
		mtf_res->conc_rand_mean = (double)conc_sum / (double)(rnd_net_num);
	else
		mtf_res->conc_rand_mean =0;

	if(rnd_net_num >1) {
		//calc counts standard deviation
		for(i=0; i<rnd_net_num; i++) {
			l_mtf=list64_get(res_tbl->rand_arr[i],(int64)mtf_res->id);
			if (l_mtf==NULL)
				x=0;
			else {
				mtf=(Motif*)l_mtf->p;
				x=mtf->count;
			}
			if( (x-mtf_res->rand_mean) == 0)
				power = 0;
			else
				power=pow((double)(x - mtf_res->rand_mean), (double)2);
			sigma += power/(double)(rnd_net_num-1);
		}
	}
	mtf_res->rand_std_dev=(double)sqrt(sigma);

	sigma=0;
	//calc concentrations standard deviation
 	for(i=0; i<rnd_net_num; i++) {
		l_mtf=list64_get(res_tbl->rand_arr[i],(int64)mtf_res->id);
		if (l_mtf==NULL)
			y=0;
		else {
			mtf=(Motif*)l_mtf->p;
			y=mtf->conc;
		}
		if( (y-mtf_res->conc_rand_mean) == 0)
			power = 0;
		else
			power=pow((double)(y - mtf_res->conc_rand_mean), (double)2);
		sigma += power/(double)(rnd_net_num-1);
	}
	if(rnd_net_num > 0)
		mtf_res->conc_rand_std_dev=(double)sqrt(sigma);
	else
		mtf_res->conc_rand_std_dev=0;

}


void
calc_pval(Motif_res *mtf_res, Res_tbl *res_tbl,int rnd_net_num)
{
	int i;
	list64_item *l_mtf;
	Motif *mtf;
	int real_pval=0;
	int real_conc_pval=0;

	//calc real pval
	for(i=0; i<rnd_net_num; i++) {
		if ( (l_mtf=list64_get(res_tbl->rand_arr[i],(int64)mtf_res->id)) !=NULL) {
			mtf=(Motif*)l_mtf->p;
			if(mtf_res->real_count<=(double)(mtf->count))
				real_pval++;
			if(mtf_res->conc_real<=(double)(mtf->conc))
				real_conc_pval++;
		//rand net has no occurances
			//if real net is zero appearances too then increase pval count by one
		} else {
			if(mtf_res->real_count==(double)0)
				real_pval++;
			if(mtf_res->conc_real==(double)0)
				real_conc_pval++;
		}

	}
	if(rnd_net_num>0){
		mtf_res->real_pval = (double)real_pval / (double)rnd_net_num;
		mtf_res->conc_real_pval = (double)real_conc_pval / (double)rnd_net_num;
	} else {
		mtf_res->real_pval=0;
		mtf_res->conc_real_pval=0;
	}
}

void
calc_zscore(Motif_res *mtf_res)
{
	if (mtf_res->rand_std_dev > 0) {
		mtf_res->real_zscore = (double)(mtf_res->real_count - mtf_res->rand_mean)
			/ (double)mtf_res->rand_std_dev;
		mtf_res->conc_real_zscore = (double)(mtf_res->conc_real - mtf_res->conc_rand_mean)
			/ (double)mtf_res->conc_rand_std_dev;
	} else {
		mtf_res->real_zscore = UNDEFINED;
		mtf_res->conc_real_zscore = UNDEFINED;
	}
}



void
calc_roots_leafs_num(Network *N)
{
	int i;

	N->roots_num=0;
	N->leafs_num=0;
	for(i=1;i<=N->vertices_num;i++) {
		if(N->indeg[i]==0 && N->outdeg[i]>0)
			N->roots_num++;
		if(N->indeg[i]>0 && N->outdeg[i]==0)
			N->leafs_num++;
	}
}

int
is_dangling(int64 id, int sz)
{
	Matrix *M,*M_T,*S;
	int i,j;
	int row_sum=0;
	int dangling=FALSE;

	M=init_matrix(sz);
	M_T=init_matrix(sz);
	S=init_matrix(sz);
	fill_mat_id(M,id);
	//M_T is the M transpose
	for(i=1;i<=sz;i++)
		for(j=1;j<=sz;j++)
			MTRX(M_T,j,i)=MTRX(M,i,j);
	//S cells are the OR of all twin cells of M and M_T
	for(i=1;i<=sz;i++)
		for(j=1;j<=sz;j++)
			MTRX(S,i,j)=(char)(MTRX(M,i,j) || MTRX(M_T,i,j));
	for(i=1;i<=sz;i++) {
		row_sum=0;
		for(j=1;j<=sz;j++)
			row_sum+=(char)MTRX(S,i,j);
		if(row_sum<2) {
			dangling=TRUE;
			break;
		}
	}
	//dump_matrix(stdout,M, "");
	//dump_matrix(stdout,M_T, "");
	//dump_matrix(stdout,S, "");
	//fprintf(stdout,"\n dangling: %d\n",dangling);

	free_matrix(M);free_matrix(M_T);free_matrix(S);
	return dangling;
}


int
calc_final_results(Res_tbl*res_tbl, list64 **final_res_p, list64 **final_res_all_p, int rnd_net_num)
{
	int rc=RC_OK;
	list64_item *l_id_res,*l_tmp;
	Motif_res *mtf_res;
	list64 *final_res, *all_ids_res;
	list64 *all_pos_ids;
	list64* iso_list;
	int64 id,rep_id,mask_bit,mask_r,mask_c;
	int i,j,illegal_id=FALSE;


	list64_init(&final_res);

	//for each motif in real network
	for(l_id_res=list64_get_next(res_tbl->real,NULL); l_id_res !=NULL;
		l_id_res=list64_get_next(res_tbl->real,l_id_res)) {
			mtf_res=(Motif_res*)calloc(1,sizeof(Motif_res));
			if(GNRL_ST.out_rand_mat_flag)
				mtf_res->all_rand_counts=(int*)calloc(GNRL_ST.rnd_net_num+1,sizeof(int));
		
			mtf_res->id=(int64)l_id_res->val;
			mtf_res->real_count=(double)((Motif*)(l_id_res->p))->count;
			mtf_res->conc_real=((Motif*)(l_id_res->p))->conc;
			mtf_res->hits_num=((Motif*)(l_id_res->p))->hits;
			mtf_res->dangling=FALSE;


			//full set motif results
			calc_deviation(mtf_res,&RES_TBL,rnd_net_num);
			calc_pval(mtf_res,&RES_TBL,rnd_net_num);
			calc_zscore(mtf_res);
			mtf_res->unique_appear=((Motif*)(l_id_res->p))->members->size;
			mtf_res->hits_num=(int)((Motif*)(l_id_res->p))->hits;
			//if list all members flag then add reference to list of all members
			if(GNRL_ST.list_members==TRUE) {
				mtf_res->all_members=((Motif*)(l_id_res->p))->all_members;
			}
			if(GNRL_ST.prob_converge_mode)
				mtf_res->conv_grade=((Motif*)(l_id_res->p))->conv_grade;

			//check if there are dangling edges in this subgraph
			mtf_res->dangling=is_dangling(mtf_res->id, GNRL_ST.mtf_sz);
			list64_insert(final_res,(int64)mtf_res->id,(void*)mtf_res);
	}
	*final_res_p=final_res;

	//if motif size is larger than 4 then skip this (takes too much time)
	if (GNRL_ST.mtf_sz<=4){
		//generate a list with representatives for all possible motifs id for
		//the motif size and then pass through the list and get all the results
		list64_init(&all_pos_ids);
		for(id=0;id<=(int64)(pow(2,GNRL_ST.mtf_sz*GNRL_ST.mtf_sz)-1);id++) {
			illegal_id=FALSE;
			//check if legal id,
			if(GNRL_ST.calc_self_edges == FALSE)
			{
				//First check if do not contain self edges
				for(i=0;i<GNRL_ST.mtf_sz;i++) {
					//bits on the diagonal
					mask_bit=(int64)pow(2,i*GNRL_ST.mtf_sz+i);
					if(id & mask_bit) {
						illegal_id = TRUE;
						break;
					}
				}
			}

			//second check there is no isolated vertex (a vertexs with no edges at all)
			//this is done by checking that there is at list one edge
			//at each row i or column i of the matrix
			//(by bit manipulation)
			for(i=0;i<GNRL_ST.mtf_sz;i++) {
				mask_r=0;
				mask_c=0;
				for(j=0;j<GNRL_ST.mtf_sz;j++) {
					mask_r |= (int64)pow(2,i*GNRL_ST.mtf_sz+j);
					mask_c |= (int64)pow(2,i+j*GNRL_ST.mtf_sz);
				}
				if( ! ((id & mask_r) || (id & mask_c)) ) {
					illegal_id = TRUE;
					break;
				}
				if(GNRL_ST.calc_self_edges == TRUE)
				{
					if(((id & mask_r) == (int64)pow(2,i*GNRL_ST.mtf_sz+i)) &&
					   ((id & mask_c) == (int64)pow(2,i+i*GNRL_ST.mtf_sz)))
					{
						illegal_id = TRUE;
						break;
					}
				}
			}
			if(illegal_id == TRUE)
				continue;
			iso_list=calc_mtf_id_iso(id,GNRL_ST.mtf_sz);
			l_tmp=list64_get_next(iso_list,NULL);
			if(l_tmp!=NULL) {
				rep_id=l_tmp->val;
				if(list64_get(all_pos_ids,(int64)rep_id)==NULL)
					//third check - single connected component
					//at the moment just remove the non-connected manually (easier)
					if( !((GNRL_ST.mtf_sz==4) && ((rep_id==72) || (rep_id==388) || (rep_id==4680))) )
						list64_insert(all_pos_ids,(int64)rep_id,NULL);
			}
			list64_free_mem(iso_list);
		}
		//pass through all possible rep ids
		//collect the result in all_ids_res
		list64_init(&all_ids_res);
		for(l_tmp=list64_get_next(all_pos_ids,NULL); l_tmp !=NULL;
		l_tmp=list64_get_next(all_pos_ids,l_tmp)) {
			rep_id=l_tmp->val;
			mtf_res=(Motif_res*)calloc(1,sizeof(Motif_res));
			if(GNRL_ST.out_rand_mat_flag)
				mtf_res->all_rand_counts=(int*)calloc(GNRL_ST.rnd_net_num+1,sizeof(int));
		
			mtf_res->id=(int64)rep_id;
			//if found in real network then get its score
			//else the score is zero
			if( (l_id_res=list64_get(res_tbl->real,(int64)rep_id)) != NULL) {
				mtf_res->real_count=(double)((Motif*)(l_id_res->p))->count;
				mtf_res->hits_num=(int)((Motif*)(l_id_res->p))->count;
				mtf_res->conc_real=((Motif*)(l_id_res->p))->conc;
				mtf_res->hits_num=((Motif*)(l_id_res->p))->hits;
				mtf_res->all_members=((Motif*)(l_id_res->p))->all_members;
				mtf_res->unique_appear=((Motif*)(l_id_res->p))->members->size;
				mtf_res->conv_grade=((Motif*)(l_id_res->p))->conv_grade;

			}else{
				mtf_res->real_count=0;
				mtf_res->conc_real=0;
			}

			//full set motif results
			calc_deviation(mtf_res,&RES_TBL,rnd_net_num);
			calc_pval(mtf_res,&RES_TBL,rnd_net_num);
			calc_zscore(mtf_res);

			list64_insert(all_ids_res,(int64)mtf_res->id,(void*)mtf_res);
		}

		*final_res_all_p=all_ids_res;
	}
	if (GNRL_ST.mtf_sz<=4)
		list64_free_mem(all_pos_ids);
	//calculate roles final results
	//currently only ofr motifs size 3
	if(GNRL_ST.calc_roles && GNRL_ST.mtf_sz==3 ) {
		calc_roles_final_res(rnd_net_num,final_res_all);
	}
	//calc roots leafs num
	calc_roots_leafs_num(G_N);
	return rc;
}

/**************Grassberger stuff***************************/

void
calc_deviation_weighted(Motif_res *mtf_res, Res_tbl *res_tbl, int rnd_net_num,double *weights_arr)
{
	int i;
	list64_item *l_mtf;
	Motif *mtf;
	double x,y;
	double sigma=0;
	double power;

	mtf_res->rand_mean=0;
	// calc counts mean val
	for(i=0; i<rnd_net_num; i++) {
		if( (l_mtf=list64_get(res_tbl->rand_arr[i],(int64)mtf_res->id)) !=NULL) {
			mtf=(Motif*)l_mtf->p;
			mtf_res->rand_mean += mtf->count*weights_arr[i];
			if(GNRL_ST.out_rand_mat_flag)
				mtf_res->all_rand_counts[i]=(int)mtf->count;
		}
	}

	// calc concentration mean val
	mtf_res->conc_rand_mean=0;
	for(i=0; i<rnd_net_num; i++) {
		if( (l_mtf=list64_get(res_tbl->rand_arr[i],(int64)mtf_res->id)) !=NULL) {
			mtf=(Motif*)l_mtf->p;
			mtf_res->conc_rand_mean += mtf->conc*weights_arr[i];
		}
	}


	if(rnd_net_num >1) {
		//calc counts standard deviation
		for(i=0; i<rnd_net_num; i++) {
			l_mtf=list64_get(res_tbl->rand_arr[i],(int64)mtf_res->id);
			if (l_mtf==NULL)
				x=0;
			else {
				mtf=(Motif*)l_mtf->p;
				x=mtf->count;
			}
			power=pow(x-mtf_res->rand_mean,(double)2);

			sigma += power*weights_arr[i];
		}
	}
	mtf_res->rand_std_dev=(double)sqrt(sigma);

	sigma=0;
	//calc concentrations standard deviation
 	for(i=0; i<rnd_net_num; i++) {
		l_mtf=list64_get(res_tbl->rand_arr[i],(int64)mtf_res->id);
		if (l_mtf==NULL)
			y=0;
		else {
			mtf=(Motif*)l_mtf->p;
			y=mtf->conc;
		}
		if( (y-mtf_res->conc_rand_mean) == 0)
			power = 0;
		else
			power=pow((double)(y - mtf_res->conc_rand_mean), (double)2);
		sigma += power*(weights_arr[i]);
	}
	if(rnd_net_num > 0)
		mtf_res->conc_rand_std_dev=(double)sqrt(sigma);
	else
		mtf_res->conc_rand_std_dev=0;

}


void
calc_pval_weighted(Motif_res *mtf_res, Res_tbl *res_tbl,int rnd_net_num,double *weights_arr)
{
	int i;
	list64_item *l_mtf;
	Motif *mtf;
	double real_pval=0;
	double real_conc_pval=0;

	//calc real pval
	//weighted version
	for(i=0; i<rnd_net_num; i++) {
		if ( (l_mtf=list64_get(res_tbl->rand_arr[i],(int64)mtf_res->id)) !=NULL) {
			mtf=(Motif*)l_mtf->p;
			if(mtf_res->real_count<=(double)(mtf->count))
				real_pval+=weights_arr[i];
			if(mtf_res->conc_real<=(double)(mtf->conc))
				real_conc_pval+=weights_arr[i];
		}
	}
}


int
calc_final_results_grassberger(Res_tbl*res_tbl, int sub_flag, list64 *res_sub_mtf, list64 **final_res_p, list64 **final_res_all_p, int rnd_net_num, double *weights_arr)
{
	int rc=RC_OK;
	list64_item *l_id_res,*l_tmp;
	Motif_res *mtf_res;
	list64 *final_res, *all_ids_res;
	list64 *all_pos_ids;
	list64* iso_list;
	int64 id,rep_id,mask_bit,mask_r,mask_c;
	int i,j,illegal_id=FALSE;


	list64_init(&final_res);

	//for each motif in real network
	for(l_id_res=list64_get_next(res_tbl->real,NULL); l_id_res !=NULL;
		l_id_res=list64_get_next(res_tbl->real,l_id_res)) {
			mtf_res=(Motif_res*)calloc(1,sizeof(Motif_res));
			if(GNRL_ST.out_rand_mat_flag)
				mtf_res->all_rand_counts=(int*)calloc(GNRL_ST.rnd_net_num+1,sizeof(int));
			
			mtf_res->id=(int64)l_id_res->val;
			mtf_res->real_count=(double)((Motif*)(l_id_res->p))->count;
			mtf_res->conc_real=((Motif*)(l_id_res->p))->conc;
			mtf_res->hits_num=((Motif*)(l_id_res->p))->hits;

			//full set motif results
			calc_deviation_weighted(mtf_res,&RES_TBL,rnd_net_num,weights_arr);
			calc_pval_weighted(mtf_res,&RES_TBL,rnd_net_num,weights_arr);
			calc_zscore(mtf_res);
			mtf_res->unique_appear=((Motif*)(l_id_res->p))->members->size;
			mtf_res->hits_num=(int)((Motif*)(l_id_res->p))->hits;
			//if list all members flag then add reference to list of all members
			if(GNRL_ST.list_members==TRUE) {
				mtf_res->all_members=((Motif*)(l_id_res->p))->all_members;
			}
			if(GNRL_ST.prob_converge_mode)
				mtf_res->conv_grade=((Motif*)(l_id_res->p))->conv_grade;


			list64_insert(final_res,(int64)mtf_res->id,(void*)mtf_res);
	}
	*final_res_p=final_res;

	//if was called by process submotifs then skip the all ids calc results
	//if motif size is greater then 4 then skip this (takes too much time)
	if ((sub_flag == FALSE) && (GNRL_ST.mtf_sz<=4)){
		//generate a list with representatives for all possible motifs id for
		//the motif size and then pass through the list and get all the results
		list64_init(&all_pos_ids);
		for(id=0;id<=(int64)(pow(2,GNRL_ST.mtf_sz*GNRL_ST.mtf_sz)-1);id++) {
			illegal_id=FALSE;
			//check if legal id,
			if(GNRL_ST.calc_self_edges == FALSE)
			{
				//First check if do not contain self edges
				for(i=0;i<GNRL_ST.mtf_sz;i++) {
					//bits on the diagonal
					mask_bit=(int64)pow(2,i*GNRL_ST.mtf_sz+i);
					if(id & mask_bit) {
						illegal_id = TRUE;
						break;
					}
				}
			}

			//second check there is no isolated vertex (a vertexs with no edges at all)
			//this is done by checking that there is at list one edge
			//at each row i or column i of the matrix
			//(by bit manipulation)
			for(i=0;i<GNRL_ST.mtf_sz;i++) {
				mask_r=0;
				mask_c=0;
				for(j=0;j<GNRL_ST.mtf_sz;j++) {
					mask_r |= (int64)pow(2,i*GNRL_ST.mtf_sz+j);
					mask_c |= (int64)pow(2,i+j*GNRL_ST.mtf_sz);
				}
				if( ! ((id & mask_r) || (id & mask_c)) ) {
					illegal_id = TRUE;
					break;
				}
				if(GNRL_ST.calc_self_edges == TRUE)
				{
					if(((id & mask_r) == (int64)pow(2,i*GNRL_ST.mtf_sz+i)) &&
					   ((id & mask_c) == (int64)pow(2,i+i*GNRL_ST.mtf_sz)))
					{
						illegal_id = TRUE;
						break;
					}
				}
			}
			//third check - connectivity
			//nadav - to add here

			if(illegal_id == TRUE)
				continue;
			iso_list=calc_mtf_id_iso(id,GNRL_ST.mtf_sz);
			l_tmp=list64_get_next(iso_list,NULL);
			if(l_tmp!=NULL) {
				rep_id=l_tmp->val;
				if(list64_get(all_pos_ids,(int64)rep_id)==NULL)
					//third check - connectivity
					//temporary :nadav just for 4 size motifs
					//- to be removed when I write a proper general third check
					if( !((GNRL_ST.mtf_sz==4) && ((rep_id==72) || (rep_id==388) || (rep_id==4680))) )
						list64_insert(all_pos_ids,(int64)rep_id,NULL);
			}
			list64_free_mem(iso_list);
			//free((void*)iso_list);
		}
		//pass through all possible rep ids
		//collect the result in all_ids_res
		list64_init(&all_ids_res);
		for(l_tmp=list64_get_next(all_pos_ids,NULL); l_tmp !=NULL;
		l_tmp=list64_get_next(all_pos_ids,l_tmp)) {
			rep_id=l_tmp->val;
			mtf_res=(Motif_res*)calloc(1,sizeof(Motif_res));
			if(GNRL_ST.out_rand_mat_flag)
				mtf_res->all_rand_counts=(int*)calloc(GNRL_ST.rnd_net_num+1,sizeof(int));
			
			mtf_res->id=(int64)rep_id;
			//if found in real network then get its score
			//else the score is zero
			if( (l_id_res=list64_get(res_tbl->real,(int64)rep_id)) != NULL) {
				mtf_res->real_count=(double)((Motif*)(l_id_res->p))->count;
				mtf_res->hits_num=(int)((Motif*)(l_id_res->p))->count;
				mtf_res->conc_real=((Motif*)(l_id_res->p))->conc;
				mtf_res->hits_num=((Motif*)(l_id_res->p))->hits;
				mtf_res->all_members=((Motif*)(l_id_res->p))->all_members;
				mtf_res->unique_appear=((Motif*)(l_id_res->p))->members->size;
				mtf_res->conv_grade=((Motif*)(l_id_res->p))->conv_grade;

			}else{
				mtf_res->real_count=0;
				mtf_res->conc_real=0;
			}

			//full set motif results
			calc_deviation_weighted(mtf_res,&RES_TBL,rnd_net_num,weights_arr);
			calc_pval_weighted(mtf_res,&RES_TBL,rnd_net_num,weights_arr);
			calc_zscore(mtf_res);

			list64_insert(all_ids_res,(int64)mtf_res->id,(void*)mtf_res);
		}

		*final_res_all_p=all_ids_res;
	}
	list64_free_mem(all_pos_ids);
	//calculate roles final results
	//currently only ofr motifs size 3
	if(GNRL_ST.calc_roles && GNRL_ST.mtf_sz==3 && !sub_flag) {
		calc_roles_final_res(rnd_net_num,final_res_all);
	}
	//calc roots leafs num
	calc_roots_leafs_num(G_N);
	return rc;
}



