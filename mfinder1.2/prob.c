/************************************************************************
*
*  File name: prob.c
*
*  Description: Probabilistic algorithm (Sampling method) functions
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
#include "motif_ids.h"
#include "metropolis.h"
#include "prob.h"
#include "results.h"





/*************************** Global variables ****************************/
char *ngbr_edges_vec;
int *ngbr_edges_indices;
char *hub_edges_vec;
int *hub_edges_indices;



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

//moved from main
extern list64 *final_res, *final_res_all;
extern list64 *res_sub_motif;



extern void init_random_seed();

// return random integer in the interval 1 to max_val
extern int get_rand(int max_val);
extern double get_rand_double();
//Roles staff
extern Hash *Role_hash;
extern role_st *Role_trans;

/******************************* Declerations ****************************/




/******************************* Functions *******************************/




//retrun degree of vertex (in + out)
int
get_deg(Network *N,int vertex)
{
	int deg=0;
	int i;
	if(N->mat->type==FULL) {
		for(i=1;i<=N->vertices_num;i++) {
			deg+=MatGet(N->mat,vertex,i);
			deg+=MatGet(N->mat,i,vertex);
		}
	}else{
		//sparse - just the sum of to and from lists
		if(N->mat->spr->m[vertex].to !=NULL)
			deg+=N->mat->spr->m[vertex].to->size;
		if(N->mat->spr->m[vertex].from !=NULL)
			deg+=N->mat->spr->m[vertex].from->size;
	}
	return deg;
}

//get effective degree of node
//this is done for the prob approach
//vertex set is the current set already chosen
int
get_eff_deg(Network *N, int vertex, list*vrtx_set)
{
	int e_deg;
	list_item *l_v;
	e_deg=get_deg(N,vertex);
	for(l_v=list_get_next(vrtx_set,NULL);l_v!=NULL;l_v=list_get_next(vrtx_set,l_v)) {
		e_deg-=MatGet(N->mat,vertex,l_v->val);
		e_deg-=MatGet(N->mat,l_v->val,vertex);
	}
	return e_deg;
}

//get effective degree of node
//this is done for the prob approach
//vertex set is the current set already chosen
int
get_eff_deg_efc(Network *N, int vertex, list*vrtx_set,Network*SN,int*sn_node_map)
{
	int e_deg;
	list_item *l_v;
	e_deg=get_deg(N,vertex);
	for(l_v=list_get_next(vrtx_set,NULL);l_v!=NULL;l_v=list_get_next(vrtx_set,l_v)) {
		e_deg-=MatGet(SN->mat,sn_node_map[vertex],sn_node_map[l_v->val]);
		e_deg-=MatGet(SN->mat,sn_node_map[l_v->val],sn_node_map[vertex]);
	}
	return e_deg;
}

//get_p_i_
//called by get_p_i in a recursive way in order to calculate p_i
double
get_p_i_(Network *N, list *vrtx_set, list *cur_set)
{
	list_item *l_v,*c_v;
	int v;
	double p_i_=0;
	double cur_p_=0;
	int nom,dom;

	//if no more vertices to complete vrtx_set
	if(cur_set->size==vrtx_set->size)
		return 1;

	for(l_v=list_get_next(vrtx_set,NULL);l_v!=NULL;l_v=list_get_next(vrtx_set,l_v)) {
		// if not a vertex in cur_set
		if(list_get(cur_set,l_v->val)==NULL) {
			v=l_v->val;

			//denominator
			//sigma e_deg for all vertices in cur_set
			dom=0;
			for(c_v=list_get_next(cur_set,NULL);c_v!=NULL;c_v=list_get_next(cur_set,c_v)) {
				dom+=get_eff_deg(N,c_v->val,cur_set);
				}
			//nominator
			//all edges from cur_set to v
			nom=0;
			for(c_v=list_get_next(cur_set,NULL);c_v!=NULL;c_v=list_get_next(cur_set,c_v)) {
				nom+=MatGet(N->mat,c_v->val,v);
				nom+=MatGet(N->mat,v,c_v->val);
			}
			if(DEBUG_LEVEL >=5 && GNRL_ST.out_log_flag)
				fprintf(GNRL_ST.log_fp,"nom = %f,dom = %f\n",(double)nom,(double)dom);
			list_insert(cur_set,v,NULL);
			cur_p_=(double)(((double)nom/(double)dom)*get_p_i_(N,vrtx_set,cur_set));
			list_delete(cur_set,v);
			p_i_+=cur_p_;
		}
	}
	return p_i_;
}


//get_p_i_
//called by get_p_i in a recursive way in order to calculate p_i
double
get_p_i_efc_(Network *N, list *vrtx_set, list *cur_set,Network *SN,int* sn_node_map)
{
	list_item *l_v,*c_v;
	int v;
	double p_i_=0;
	double cur_p_=0;
	int nom,dom;

	//if no more vertices to complete vrtx_set
	if(cur_set->size==vrtx_set->size)
		return 1;

	for(l_v=list_get_next(vrtx_set,NULL);l_v!=NULL;l_v=list_get_next(vrtx_set,l_v)) {
		// if not a vertex in cur_set
		if(list_get(cur_set,l_v->val)==NULL) {
			v=l_v->val;

			//denominator
			//sigma e_deg for all vertices in cur_set
			dom=0;
			for(c_v=list_get_next(cur_set,NULL);c_v!=NULL;c_v=list_get_next(cur_set,c_v)) {
				dom+=get_eff_deg_efc(N,c_v->val,cur_set,SN,sn_node_map);
				}
			//nominator
			//all edges from cur_set to v
			nom=0;
			for(c_v=list_get_next(cur_set,NULL);c_v!=NULL;c_v=list_get_next(cur_set,c_v)) {
				nom+=MatGet(SN->mat,sn_node_map[c_v->val],sn_node_map[v]);
				nom+=MatGet(SN->mat,sn_node_map[v],sn_node_map[c_v->val]);
			}
			if(DEBUG_LEVEL >=5 && GNRL_ST.out_log_flag)
				fprintf(GNRL_ST.log_fp,"nom = %f,dom = %f\n",(double)nom,(double)dom);
			list_insert(cur_set,v,NULL);
			cur_p_=(double)(((double)nom/(double)dom)*get_p_i_efc_(N,vrtx_set,cur_set,SN,sn_node_map));
			list_delete(cur_set,v);
			p_i_+=cur_p_;
		}
	}
	return p_i_;
}

//get pi - the probablilty to reach this vertex_set
//in the probabilistic alg
double
get_p_i(Network *N, list *vrtx_set, Network *SN)
{
  //int rc=RC_OK;
	Edge e;
	list_item *l_v;
	list *cur_set;
	int i,j;
	double p_i=0;
	double cur_p=0;


	//if motif size 2 then p_i is 1/#edges
	if(vrtx_set->size==2) {
		p_i = (double)1/(double)N->edges_num;
		return p_i;
	}
	//if motif size >=3

	for(i=1;i<=vrtx_set->size;i++) {
		for(j=1;j<=vrtx_set->size;j++) {
			//if there is an edge e(i,j)
			if(MatGet(SN->mat,i,j) != 0) {

				l_v=list_get_by_indx(vrtx_set,i);
				if( (l_v==NULL) )
					return RC_ERR;
				else
					e.s=l_v->val;
				l_v=list_get_by_indx(vrtx_set,j);
				if( (l_v==NULL) )
					return RC_ERR;
				else
					e.t=l_v->val;

				list_init(&cur_set);
				list_insert(cur_set,e.s,NULL);
				list_insert(cur_set,e.t,NULL);
				if(DEBUG_LEVEL >=5 && GNRL_ST.out_log_flag)
					fprintf(GNRL_ST.log_fp,"s = %d,t= %d\n",e.s,e.t);
				cur_p=(double)((double)1/(double)N->edges_num)*get_p_i_(N,vrtx_set,cur_set);
				list_free_mem(cur_set);

				p_i+=cur_p;
			}
		}
	}
	return p_i;
}


//get pi - the probablilty to reach this vertex_set
//in the probabilistic alg
double
get_p_i_efc(Network *N, list *vrtx_set, Network *SN)
{
  //int rc=RC_OK;
	Edge e;
	list_item *l_v;
	list *cur_set;
	int i,j;
	double p_i=0;
	double cur_p=0;
	int *sn_node_map;

	//if motif size 2 then p_i is 1/#edges
	if(vrtx_set->size==2) {
		p_i = (double)1/(double)N->edges_num;
		return p_i;
	}
	//if motif size >=3

	//generate node to index mapping (the index of each node in SN)
	sn_node_map=(int*)calloc(N->vertices_num+1,sizeof(int));
	for(i=1,l_v=list_get_next(vrtx_set,NULL);l_v!=NULL;i++,l_v=list_get_next(vrtx_set,l_v))
		sn_node_map[l_v->val]=i;

	for(i=1;i<=vrtx_set->size;i++) {
		for(j=1;j<=vrtx_set->size;j++) {
			//if there is an edge e(i,j)
			if(MatGet(SN->mat,i,j) != 0) {

				l_v=list_get_by_indx(vrtx_set,i);
				if( (l_v==NULL) )
					return RC_ERR;
				else
					e.s=l_v->val;
				l_v=list_get_by_indx(vrtx_set,j);
				if( (l_v==NULL) )
					return RC_ERR;
				else
					e.t=l_v->val;

				list_init(&cur_set);
				list_insert(cur_set,e.s,NULL);
				list_insert(cur_set,e.t,NULL);
				if(DEBUG_LEVEL >=5 && GNRL_ST.out_log_flag)
					fprintf(GNRL_ST.log_fp,"s = %d,t= %d\n",e.s,e.t);
				cur_p=(double)((double)1/(double)N->edges_num)*get_p_i_efc_(N,vrtx_set,cur_set,SN,sn_node_map);
				list_free_mem(cur_set);

				p_i+=cur_p;
			}
		}
	}
	free(sn_node_map);
	return p_i;
}



//get pi - the probablilty to reach this vertex_set
//in the probabilistic alg
double
get_approx_p_i(Network *N, list *vrtx_set, Network *SN)
{
  //int rc=RC_OK;
	//Edge e;
	list_item *l_v;
	int i,j;
	double p_i=0;
	int internal_edge_num=0;
	//int internal_deg=0;
	int sigma_internal_deg=0;
	int sigma_external_deg=0;

	//if motif size 2 then p_i is 1/#edges
	if(vrtx_set->size==2) {
		p_i = (double)1/(double)N->edges_num;
		return p_i;
	}
	//if motif size >=3
	//dump_network(stdout,SN);
	for(i=1;i<=vrtx_set->size;i++) {
		for(j=1;j<=vrtx_set->size;j++) {
			//if there is an edge e(i,j)
			if(MatGet(SN->mat,i,j) != 0)
				internal_edge_num++;
		}
	}
	//each edge was counted twice
	//each edge donate 2 to the internal degree
	sigma_internal_deg=2*internal_edge_num;


	for(l_v=list_get_next(vrtx_set,NULL);l_v!=NULL;
	l_v=list_get_next(vrtx_set,l_v)){
		sigma_external_deg+=N->indeg[l_v->val];
		sigma_external_deg+=N->outdeg[l_v->val];
	}
	//P= ( sigma(internal deg of all V's) / sigma (deg of all V's) )  *  m/E
	p_i=((double)sigma_internal_deg/(double)sigma_external_deg) * ((double)internal_edge_num/(double)N->edges_num);

	return p_i;
}





/********************************************************
* function : search_subset_prob_efc
*   This is the current used version of search_subset_prob
*   (more effiecient than before)
*   Recursive function, search for subsets connected to the current vertex set
*   if vertex set equals motif size then end of recursion
*	results are returned in res list
* arguments:
*	N: Network to be searched
*	vrtx_set: list of vertices already participate in this search tree
*   ngbr_edges : list of neigbored edges
*   mtf_sz: motif size
*	res : list pointer of results
*	net_type : REAL_NET - if real network
*              RAND_NET - if random network
* return values:
*	RC_OK - if not error occured
*   RC_ERR - if error occured
*********************************************************/
int
search_subset_prob_efc(Network *N, list *vrtx_set, int mtf_sz, list64 *res, int net_type)
{
	int i,k,rc=RC_OK;
	list_item *l_node,*l_edge,*item;
	list64_item *l_res_id;
	Network *SN;
	Motif *mtf;
	int64 mtf_id,rep_mtf_id=0;
	int role_id;
	int *mtf_vrtx_arr;  //array that store the motif vertices
	Member *mtf_members_arr;
	int unique;
	int new_members;
	role_st *role;
	Edge *edge;
	int ngbr_node=0;
	double p_i;
	list_item *l_e1,*l_e2;
	int *ngbr_vec,*ngbr_indices;
	int ngbr_num,e_idx;

	//number of vertices in set is equal to motif size
	//this is the recursive base
	if (mtf_sz == vrtx_set->size) {
		//printf("found motif\n");
		gen_subset_matrix(N, &SN, vrtx_set, mtf_sz);
		//dump_network(SN);
		//calculate motif id
		mtf_id = get_sn_motif_id(SN);
		//printf("motif id %d\n", (int)mtf_id);
		//update res tbl
		//at the moment list "val" field is int means 32 bit
		//therefore limited to motif size of 5
		//will change val field to int64 - then limited to
		//motif size of 8
		if ( (l_res_id=list64_get(res,(int64)mtf_id)) == NULL ){
			mtf=(Motif*)calloc(1,sizeof(Motif));
			mtf->id = mtf_id;
			mtf->hits=1;
			p_i=get_p_i_efc(N,vrtx_set,SN);
			//p_i=get_p_i(N,vrtx_set,SN);
			if(DEBUG_LEVEL >=5 && GNRL_ST.out_log_flag)
				fprintf(GNRL_ST.log_fp,"p_i is :%f\n",p_i);
			mtf->prob_count=(double)(1/p_i);
			list_init(&mtf->members);
			list_init(&mtf->all_members);
			if( (GNRL_ST.calc_unique_flag == TRUE) && (net_type == REAL_NET) && (mtf_sz==GNRL_ST.mtf_sz)) {
				//alloc mtf_vrtx_arr copy to it the motif vertices from vrtx_set
				mtf_vrtx_arr = (int*)calloc(mtf_sz,sizeof(int));
				
				for(i=0,item=list_get_next(vrtx_set,NULL);item!=NULL;
				i++,item=list_get_next(vrtx_set,item)) {
					mtf_vrtx_arr[i]=item->val;
				}
				list_insert(mtf->members,1,(void*)mtf_vrtx_arr);
							}
			if( ((GNRL_ST.calc_roles) || ((GNRL_ST.list_members)&&(net_type == REAL_NET)) )&& (mtf_sz==GNRL_ST.mtf_sz)) {
				//alloc mtf_members_arr copy to it the motif vertices from vrtx_set
				mtf_members_arr = (Member*)calloc(mtf_sz,sizeof(Member));
				
				//get role rules for that id
				role = hash_get(Role_hash,(int)mtf_id);
				if(GNRL_ST.mtf_sz==3)
					rep_mtf_id=(int)get_rep_mtf_id_3(mtf_id);
				for(i=0,item=list_get_next(vrtx_set,NULL);item!=NULL;
				i++,item=list_get_next(vrtx_set,item)) {
					mtf_members_arr[i].node=item->val;
					//currently support roles only for size 3
					if(GNRL_ST.mtf_sz==3) {
						mtf_members_arr[i].role=role->roles[i];
						//only if calc roles
						if(GNRL_ST.calc_roles) {
							role_id=Role_trans[rep_mtf_id].roles[mtf_members_arr[i].role];
							N->roles_vec[mtf_members_arr[i].node][role_id]=1;
							if (DEBUG_LEVEL ==-10 && GNRL_ST.out_log_flag)
								fprintf(GNRL_ST.log_fp,"added to net->role_vec: role_id: %d node %d\n",role_id,mtf_members_arr[i].node);

						}
					}
					else
						mtf_members_arr[i].role=1;

				}
				if((GNRL_ST.list_members == TRUE) && (net_type == REAL_NET) ) {
					list_insert(mtf->all_members,1,(void*)mtf_members_arr);
					
				} else {
					free(mtf_members_arr);
				}
			}
			//now after mtf is fixed insert it to res list (in the p field)
			//now after mtf is fixed insert it to res list (in the p field)
			list64_insert(res,(int64)mtf_id,(void*)mtf);
		}else {
			mtf=((Motif*)(l_res_id->p));
			//increase count by one
			mtf->hits+=1;

			p_i=get_p_i_efc(N,vrtx_set,SN);
			//p_i=get_p_i(N,vrtx_set,SN);
			
			mtf->prob_count+=(double)(1/p_i);

			if( (GNRL_ST.calc_unique_flag == TRUE) && (net_type == REAL_NET) && (mtf_sz==GNRL_ST.mtf_sz)) {
				//alloc mtf_vrtx_arr copy to it the motif vertices from vrtx_set
				mtf_vrtx_arr = (int*)calloc(mtf_sz,sizeof(int));
				

				for(i=0,item=list_get_next(vrtx_set,NULL);item!=NULL;
				i++,item=list_get_next(vrtx_set,item)) {
					mtf_vrtx_arr[i]=item->val;
				}
				mtf=(Motif*)l_res_id->p;

				
				//pass through list , if this vertices group is unique (no common vertices with
				//same motif id that already found) then insert to list
				unique = check_if_unique(mtf,mtf_vrtx_arr,mtf_sz);

				if(unique==TRUE) {
					if (DEBUG_LEVEL >= 1 && GNRL_ST.out_log_flag) {
						fprintf(GNRL_ST.log_fp,"Found unique motif members, going to insert them to list\n");
					}
					list_insert(mtf->members,(int)((Motif*)(l_res_id->p))->hits,(void*)mtf_vrtx_arr);
				} else {
					free(mtf_vrtx_arr);
				}
			}
			if( ((GNRL_ST.calc_roles) || ((GNRL_ST.list_members)&&(net_type == REAL_NET)) ) && (mtf_sz==GNRL_ST.mtf_sz) ) {
				//alloc mtf_members_arr copy to it the motif vertices from vrtx_set
				mtf_members_arr = (Member*)calloc(mtf_sz,sizeof(Member));
			
				//get role rules for that id
				role = hash_get(Role_hash,(int)mtf_id);
				if(GNRL_ST.mtf_sz==3)
					rep_mtf_id=(int)get_rep_mtf_id_3(mtf_id);
				for(i=0,item=list_get_next(vrtx_set,NULL);item!=NULL;
					i++,item=list_get_next(vrtx_set,item)) {
					mtf_members_arr[i].node=item->val;
					//currently support roles only for size 3
					if(GNRL_ST.mtf_sz==3) {
						mtf_members_arr[i].role=role->roles[i];
						//only if calc roles
						if(GNRL_ST.calc_roles) {
							role_id=Role_trans[rep_mtf_id].roles[mtf_members_arr[i].role];
							N->roles_vec[mtf_members_arr[i].node][role_id]=1;
							if (DEBUG_LEVEL ==-10 && GNRL_ST.out_log_flag)
								fprintf(GNRL_ST.log_fp,"added to net->role_vec: role_id: %d node %d\n",role_id,mtf_members_arr[i].node);
						}
					}
					else
						mtf_members_arr[i].role=1;
				}
				mtf=(Motif*)l_res_id->p;
				if((GNRL_ST.list_members == TRUE) && (net_type == REAL_NET) ) {
					
					new_members = check_if_new_members(mtf,mtf_members_arr,mtf_sz);
					if(new_members==TRUE)
						list_insert(mtf->all_members, 1,(void*)mtf_members_arr);
					else
						free(mtf_members_arr);
				} else {
					free(mtf_members_arr);
				}
			}
		}
		//free subset matrix
		free_subset_matrix_mem(SN);
		free((void*)SN);
		return RC_OK;
	}

	//find neighbors off all vertices in vrtx_set
	//insert all edges to ngbr_edges list
	//then choose randomly an edge from this list

	//allocate ngbr vec and ngbr indices arr
	ngbr_vec=(int*)calloc(N->edges_num+1,sizeof(int));
	ngbr_indices=(int*)calloc((N->hub_deg+1)*GNRL_ST.mtf_sz,sizeof(int));
	ngbr_num=0;

	if(DEBUG_LEVEL ==-100 && GNRL_ST.out_log_flag) {
		fprintf(GNRL_ST.log_fp,"vertex set (size %d) :\t",vrtx_set->size);
		for(l_node=list_get_next(vrtx_set,NULL);l_node!=NULL;
			l_node=list_get_next(vrtx_set,l_node)) {
			fprintf(GNRL_ST.log_fp,"%d ",l_node->val);
		}
		fprintf(GNRL_ST.log_fp,"\n");
		fflush(GNRL_ST.log_fp);
	}
	//Check if the hub is in vrtx set if yes then take its vec and indices
	//from the already made Network structure
	for(l_node=list_get_next(vrtx_set,NULL);l_node!=NULL;
	l_node=list_get_next(vrtx_set,l_node)) {
		//if the HUB
		if(l_node->val==N->hub) {
			memcpy(ngbr_vec,N->hub_edges_vec,(N->edges_num+1)*sizeof(int));
			memcpy(ngbr_indices,N->hub_edges_indices,(N->hub_deg+1)*sizeof(int));
			ngbr_num+=N->hub_deg;

			
		}
	}

	for(l_node=list_get_next(vrtx_set,NULL);l_node!=NULL;
	l_node=list_get_next(vrtx_set,l_node)) {
		//if not the HUB
		if(l_node->val != N->hub) {
			//pass throgh e_mat, it is a sparse matrix
			//pass through row
			for(l_edge=list_get_next(N->e_map->spr->m[l_node->val].to,NULL);l_edge!=NULL;
			l_edge=list_get_next(N->e_map->spr->m[l_node->val].to,l_edge)) {
				//if a "to" neighbor of the current vertices set
				//and not yet in the neigbours edges
				e_idx=*(int*)l_edge->p;
				if(ngbr_vec[e_idx]==0) {
					//the indice of the new neighbor point to the edge index
					ngbr_indices[++ngbr_num]=e_idx;
					//the edge entry in the vec hold its neigbor number
					ngbr_vec[e_idx]=ngbr_num;
				}
			}
			//pass through column
			for(l_edge=list_get_next(N->e_map->spr->m[l_node->val].from,NULL);l_edge!=NULL;
			l_edge=list_get_next(N->e_map->spr->m[l_node->val].from,l_edge)) {
				//if a "from" neighbor of the current vertices set
				//and not yet in the neigbours edges
				e_idx=*(int*)l_edge->p;
				if(ngbr_vec[e_idx]==0) {
					//the indice of the new neighbor point to the edge index
					ngbr_indices[++ngbr_num]=e_idx;
					//the edge entry in the vec hold its neigbor number
					ngbr_vec[e_idx]=ngbr_num;
				}
			}
		}
	}//for loop on vrtx_set

	if(DEBUG_LEVEL ==-100 && GNRL_ST.out_log_flag) {
		fprintf(GNRL_ST.log_fp,"ngbr arrays BEFORE intern edges remove vertex set\n");
		fprintf(GNRL_ST.log_fp,"Ngbr_vec : ");
		for(i=1;i<=N->edges_num;i++)
			fprintf(GNRL_ST.log_fp,"%d ",ngbr_vec[i]);
		fprintf(GNRL_ST.log_fp,"\nNgbr_indices : ");
		for(i=1;i<=N->hub_deg;i++)
			fprintf(GNRL_ST.log_fp,"%d ",ngbr_indices[i]);
		fprintf(GNRL_ST.log_fp,"\n");
		fflush(GNRL_ST.log_fp);
	}

	//Now left to remove all edges that are between nodes in the current set of nodes
	//This is done very efficiently using the two arrays , a little bit fumbling
	for(l_e1=list_get_next(vrtx_set,NULL);l_e1!=NULL;
	l_e1=list_get_next(vrtx_set,l_e1)) {
		for(l_e2=list_get_next(vrtx_set,NULL);l_e2!=NULL;
		l_e2=list_get_next(vrtx_set,l_e2)) {
			if( (e_idx=MatGet(N->e_map,l_e1->val,l_e2->val)) != 0 ) {
				//remove this edge from the list
				//this operation based on two steps:
				//zero the edge entry in the ngbr vec
				//and remove indice from ngbr_indices
				//by copying the last entry to the now just deleted entry
				//and update ngbr_num--
				ngbr_indices[ngbr_vec[e_idx]]=ngbr_indices[ngbr_num];
				ngbr_vec[ngbr_indices[ngbr_num]]=ngbr_vec[e_idx];
				ngbr_indices[ngbr_num]=0;
				ngbr_vec[e_idx]=0;
				ngbr_num--;
			}
		}
	}
	

	//if no neigbours then we didnt find motif - return
	//without updating any motif counter
	if(ngbr_num == 0) {
		free(ngbr_vec);
		free(ngbr_indices);
		//printf("didnt reach a motif\t");
		return RC_ERR;
	}

	//get a random edge from ngbr_edges and recursively call search_subset_prob
	k=get_rand(ngbr_num);
	//check who is the new neighbor and add it to vrtx_set
	e_idx=ngbr_indices[k];
	edge=(Edge*)&N->e_arr[e_idx];

	free(ngbr_vec);
	free(ngbr_indices);

    
	if( (list_get(vrtx_set,edge->s)==NULL) && (list_get(vrtx_set,edge->t)!=NULL))
		ngbr_node=edge->s;
	else if ( (list_get(vrtx_set,edge->t)==NULL) && (list_get(vrtx_set,edge->s)!=NULL))
		ngbr_node=edge->t;
	else
		printf ("Error: no leg in leg of this edge\n");

	//insert the additional node to the current set of nodes
	list_insert(vrtx_set,ngbr_node,NULL);

	rc = search_subset_prob_efc(N, vrtx_set, mtf_sz, res, net_type);

	list_delete(vrtx_set,ngbr_node);

	return rc;
}




/********************************************************
* function : search_subset_prob
*   This is the original function, search_subset_prob_efc replaces it now
*   Recursive function, search for subsets connected to the current vewrtex set
*   if vertex set equals motif size then end of recursion
*	results are returned in res list
* arguments:
*	N: Network to be searched
*	vrtx_set: list of vertices already participate in this search tree
*   ngbr_edges : list of neigbored edges
*   mtf_sz: motif size
*	res : list pointer of results
*	net_type : REAL_NET - if real network
*              RAND_NET - if random network
* return values:
*	RC_OK - if not error occured
*   RC_ERR - if error occured
*********************************************************/
int
search_subset_prob(Network *N, list *vrtx_set, int mtf_sz, list64 *res, int net_type)
{
	int i,k,val,rc=RC_OK;
	list *ngbr_edges;
	list_item *l_node,*l_edge,*item;
	list64_item *l_res_id;
	Network *SN;
	Motif *mtf;
	int64 mtf_id,rep_mtf_id=0;
	int role_id;
	int *mtf_vrtx_arr;  //array that store the motif vertices
	Member *mtf_members_arr=NULL;
	int unique;
	int new_members;
	role_st *role;
	Edge *edge;
	int ngbr_node=0;
	double p_i;

	//number of vertices in set is equal to motif size
	//this is the recursive base
	if (mtf_sz == vrtx_set->size) {

		gen_subset_matrix(N, &SN, vrtx_set, mtf_sz);

		//calculate motif id
		mtf_id = get_sn_motif_id(SN);
		if ( (l_res_id=list64_get(res,(int64)mtf_id)) == NULL ){
			mtf=(Motif*)calloc(1,sizeof(Motif));
			mtf->id = mtf_id;
			mtf->hits=1;
			p_i=get_p_i(N,vrtx_set,SN);
			if(DEBUG_LEVEL >=5 && GNRL_ST.out_log_flag)
				fprintf(GNRL_ST.log_fp,"p_i is :%f\n",p_i);
			mtf->prob_count=(double)(1/p_i);
			list_init(&mtf->members);
			list_init(&mtf->all_members);
			if( (GNRL_ST.calc_unique_flag == TRUE) && (net_type == REAL_NET) && (mtf_sz==GNRL_ST.mtf_sz)) {
				//alloc mtf_vrtx_arr copy to it the motif vertices from vrtx_set
				mtf_vrtx_arr = (int*)calloc(mtf_sz,sizeof(int));
				
				for(i=0,item=list_get_next(vrtx_set,NULL);item!=NULL;
				i++,item=list_get_next(vrtx_set,item)) {
					mtf_vrtx_arr[i]=item->val;
				}
				list_insert(mtf->members,1,(void*)mtf_vrtx_arr);
				
			}
			if( ((GNRL_ST.calc_roles) || ((GNRL_ST.list_members)&&(net_type == REAL_NET)) )&& (mtf_sz==GNRL_ST.mtf_sz)) {
				//alloc mtf_members_arr copy to it the motif vertices from vrtx_set
				mtf_members_arr = (Member*)calloc(mtf_sz,sizeof(Member));
				//get role rules for that id
				role = hash_get(Role_hash,(int)mtf_id);
				if(GNRL_ST.mtf_sz==3)
					rep_mtf_id=(int)get_rep_mtf_id_3(mtf_id);
				for(i=0,item=list_get_next(vrtx_set,NULL);item!=NULL;
				i++,item=list_get_next(vrtx_set,item)) {
					mtf_members_arr[i].node=item->val;
					//currently support roles only for size 3
					if(GNRL_ST.mtf_sz==3) {
						mtf_members_arr[i].role=role->roles[i];
						//only if calc roles
						if(GNRL_ST.calc_roles) {
							role_id=Role_trans[rep_mtf_id].roles[mtf_members_arr[i].role];
							N->roles_vec[mtf_members_arr[i].node][role_id]=1;
							if (DEBUG_LEVEL ==-10 && GNRL_ST.out_log_flag)
								fprintf(GNRL_ST.log_fp,"added to net->role_vec: role_id: %d node %d\n",role_id,mtf_members_arr[i].node);

						}
					}
					else
						mtf_members_arr[i].role=1;

				}
				if((GNRL_ST.list_members == TRUE) && (net_type == REAL_NET) ) {
					list_insert(mtf->all_members,1,(void*)mtf_members_arr);
					
				} else {
					free(mtf_members_arr);
				}
			}
			//now after mtf is fixed insert it to res list (in the p field)
			list64_insert(res,(int64)mtf_id,(void*)mtf);
		}else {
			mtf=((Motif*)(l_res_id->p));
			//increase count by one
			mtf->hits+=1;
			p_i=get_p_i(N,vrtx_set,SN);
			if(DEBUG_LEVEL >=5 && GNRL_ST.out_log_flag)
				fprintf(GNRL_ST.log_fp,"p_i is :%f\n",p_i);
			mtf->prob_count+=(double)(1/p_i);

			if( (GNRL_ST.calc_unique_flag == TRUE) && (net_type == REAL_NET) && (mtf_sz==GNRL_ST.mtf_sz)) {
				//alloc mtf_vrtx_arr copy to it the motif vertices from vrtx_set
				mtf_vrtx_arr = (int*)calloc(mtf_sz,sizeof(int));
				

				for(i=0,item=list_get_next(vrtx_set,NULL);item!=NULL;
				i++,item=list_get_next(vrtx_set,item)) {
					mtf_vrtx_arr[i]=item->val;
				}
				mtf=(Motif*)l_res_id->p;

				
				//pass through list , if this vertices group is unique (no common vertices with
				//same motif id that already found) then insert to list
				unique = check_if_unique(mtf,mtf_vrtx_arr,mtf_sz);

				if(unique==TRUE) {
					if (DEBUG_LEVEL >= 1 && GNRL_ST.out_log_flag) {
						fprintf(GNRL_ST.log_fp,"Found unique motif members, going to insert them to list\n");
					}
					list_insert(mtf->members,(int)((Motif*)(l_res_id->p))->hits,(void*)mtf_vrtx_arr);
				} else {
					free(mtf_vrtx_arr);
				}
			}
			if( ((GNRL_ST.calc_roles) || ((GNRL_ST.list_members)&&(net_type == REAL_NET)) ) && (mtf_sz==GNRL_ST.mtf_sz) ) {
				//alloc mtf_members_arr copy to it the motif vertices from vrtx_set
				mtf_members_arr = (Member*)calloc(mtf_sz,sizeof(Member));
				
				//get role rules for that id
				role = hash_get(Role_hash,(int)mtf_id);
				if(GNRL_ST.mtf_sz==3)
					rep_mtf_id=(int)get_rep_mtf_id_3(mtf_id);
				for(i=0,item=list_get_next(vrtx_set,NULL);item!=NULL;
					i++,item=list_get_next(vrtx_set,item)) {
					mtf_members_arr[i].node=item->val;
					//currently support roles only for size 3
					if(GNRL_ST.mtf_sz==3) {
						mtf_members_arr[i].role=role->roles[i];
						//only if calc roles
						if(GNRL_ST.calc_roles) {
							role_id=Role_trans[rep_mtf_id].roles[mtf_members_arr[i].role];
							N->roles_vec[mtf_members_arr[i].node][role_id]=1;
							if (DEBUG_LEVEL ==-10 && GNRL_ST.out_log_flag)
								fprintf(GNRL_ST.log_fp,"added to net->role_vec: role_id: %d node %d\n",role_id,mtf_members_arr[i].node);
						}
					}
					else
						mtf_members_arr[i].role=1;
				}
				mtf=(Motif*)l_res_id->p;
				if((GNRL_ST.list_members == TRUE) && (net_type == REAL_NET) ) {
					
					if (DEBUG_LEVEL >= 1 && GNRL_ST.out_log_flag) {
						fprintf(GNRL_ST.log_fp,"Found unique motif members, going to insert them to list\n");
					}
					new_members = check_if_new_members(mtf,mtf_members_arr,mtf_sz);
					if(new_members==TRUE)
						list_insert(mtf->all_members, 1,(void*)mtf_members_arr);
					else
						free(mtf_members_arr);
				} else {
					free(mtf_members_arr);
				}
			}
		}
		//free subset matrix
		free_subset_matrix_mem(SN);
		free((void*)SN);
		return RC_OK;
	}

	//find neighbors off all vertices in vrtx_set
	//insert all edges to ngbr_edges list
	//then choose randomly an edge from this list
	list_init(&ngbr_edges);
	for(l_node=list_get_next(vrtx_set,NULL);l_node!=NULL;
	l_node=list_get_next(vrtx_set,l_node)) {

		//if full matrix then pass through row and column
		if(N->mat->type == FULL) {
			for(i=1;i<=N->vertices_num;i++) {
				//if a neighbor of the current vertices set
				if(MatGet(N->mat,l_node->val,i)==1) {
					//and not in vrtx_set already
					if(list_get(vrtx_set,i)==NULL) {
						//insert to ngbr_edges
						edge=(Edge*)calloc(1,sizeof(Edge));
						edge->s=l_node->val;
						edge->t=i;
						edge->weight=1;

						list_insert(ngbr_edges,ngbr_edges->size+1,(void*)edge);
					}
				}
			}
			//pass through column
			for(i=1;i<=N->vertices_num;i++) {
				//if a neighbor of the current vertices set
				if(MatGet(N->mat,i,l_node->val)==1) {
					//and not in vrtx_set already
					if(list_get(vrtx_set,i)==NULL) {
						//insert to ngbr_edges
						edge=(Edge*)calloc(1,sizeof(Edge));
						edge->s=i;
						edge->t=l_node->val;
						edge->weight=1;

						list_insert(ngbr_edges,ngbr_edges->size+1,(void*)edge);
					}
				}
			}
		} else {
			//sparse matrix
			//pass through row
			for(l_edge=list_get_next(N->mat->spr->m[l_node->val].to,NULL);l_edge!=NULL;
			l_edge=list_get_next(N->mat->spr->m[l_node->val].to,l_edge)) {
				//if a "to" neighbor of the current vertices set
				//and not in vrtx_set already
				if(list_get(vrtx_set,l_edge->val)==NULL) {
					//insert to ngbr_edges
					edge=(Edge*)calloc(1,sizeof(Edge));
					edge->s=l_node->val;
					edge->t=l_edge->val;
					edge->weight=1;

					list_insert(ngbr_edges,ngbr_edges->size+1,(void*)edge);
				}
			}
			//pass through column
			for(l_edge=list_get_next(N->mat->spr->m[l_node->val].from,NULL);l_edge!=NULL;
			l_edge=list_get_next(N->mat->spr->m[l_node->val].from,l_edge)) {
				//if a "from" neighbor of the current vertices set
				//and not in vrtx_set already
				if(list_get(vrtx_set,l_edge->val)==NULL) {
					//insert to ngbr_edges
					edge=(Edge*)calloc(1,sizeof(Edge));
					edge->s=l_edge->val;
					edge->t=l_node->val;
					edge->weight=1;
					list_insert(ngbr_edges,ngbr_edges->size+1,(void*)edge);
				}
			}
		} //sparse matrix
	}//for loop on vrtx_set
	

	//if no neigbours then we didnt find motif - return
	//without updating any motif counter
	if(ngbr_edges->size == 0) {
		//printf("didnt reach a motif\t");
		return RC_ERR;
	}

	//remove from ngbr_edges all edges that now have both source and target in vrtx_set
	for(l_edge=list_get_next(ngbr_edges,NULL);l_edge!=NULL;) {
		edge=(Edge*)l_edge->p;
		if( (list_get(vrtx_set,edge->s)!=NULL) && (list_get(vrtx_set,edge->t)!=NULL) ) {
			val=l_edge->val;
			l_edge=l_edge->next;
			list_delete(ngbr_edges,val);
		} else {
			l_edge=l_edge->next;
		}
	}

	//get a random edge from ngbr_edges and recursively call search_subset_prob
	k=get_rand(ngbr_edges->size);
	//check who is the new neighbor and add it to vrtx_set
	l_edge=list_get_by_indx(ngbr_edges,k);
	edge=(Edge*)l_edge->p;
	
	if( (list_get(vrtx_set,edge->s)==NULL) && (list_get(vrtx_set,edge->t)!=NULL))
		ngbr_node=edge->s;
	else if ( (list_get(vrtx_set,edge->t)==NULL) && (list_get(vrtx_set,edge->s)!=NULL))
		ngbr_node=edge->t;
	else
		printf ("Error: no leg in leg of this edge\n");
	list_free_mem(ngbr_edges);
	list_insert(vrtx_set,ngbr_node,NULL);
	if(GNRL_ST.efc_prob)
		rc = search_subset_prob_efc(N, vrtx_set, mtf_sz, res, net_type);
	else
		rc = search_subset_prob(N, vrtx_set, mtf_sz, res, net_type);

	list_delete(vrtx_set,ngbr_node);

	return rc;
}




void
build_new_conc_list(list64 **new_conc_list_p,list64 *res)
{
	Motif *mtf;
	list64_item *res_item;
	list64 *new_conc_list;
	double *conc_p;

	list64_init(&new_conc_list);

	if(res != NULL) {
		for(res_item=list64_get_next(res,NULL);res_item!=NULL;res_item=list64_get_next(res,res_item)) {
			mtf = (Motif*)res_item->p;
			conc_p=(double*)calloc(1,sizeof(double));
			*conc_p=mtf->conc;
			//insert to new_conc_list
			//val is the id (int64)
			//the p entry points to the conc (type :double)
			list64_insert(new_conc_list,res_item->val,(void*)conc_p);
		}
	}
	*new_conc_list_p=new_conc_list;
}

void
calc_diff_conc_lists(list64 *new_conc_list,list64 *old_conc_list,double *avg_diff_p,double *max_diff_p,list64 *res)
{
	double diff=0,avg_diff=0,max_diff=0,total_diff=0;
	list64_item *new_item,*old_item,*l_i;
	double new_conc,old_conc;
	int over_th_subgraphs_num=0;

	//path throught the lists together

	//the traversing is based oin the fact that
	//new list contains all subgraphs ids which the old contains
	for(new_item=list64_get_next(new_conc_list,NULL);new_item!=NULL;
	new_item=list64_get_next(new_conc_list,new_item)) {
		new_conc=*(double*)new_item->p;
		//if subthreshold conc then do not take into acount in the
		//diff calculations
		if(new_conc>=GNRL_ST.prob_conv_conc_th)
			over_th_subgraphs_num++;
		old_item=list64_get(old_conc_list,new_item->val);
		if(old_item==NULL)
			old_conc=0;
		else
			old_conc=*(double*)old_item->p;
		//update diff (add absolut value)
		if ( (new_conc-old_conc) >= 0 )
			diff=(new_conc-old_conc)/((new_conc+old_conc)/2);
		else
			diff=(old_conc-new_conc)/((new_conc+old_conc)/2);
		if(new_conc>=GNRL_ST.prob_conv_conc_th) {
			//update avg diff and max diff
			total_diff+=diff;
			if (diff>max_diff)
				max_diff=diff;
		}
		//insert conv diff to res list
		if( (l_i=list64_get(res,new_item->val))!=NULL)
			((Motif*)(l_i->p))->conv_grade=diff;
	}
	//left to normalise avg diff by the number of subgraphs
	avg_diff=total_diff/(double)over_th_subgraphs_num;
	fprintf(stdout,"Subgraphs sampled: Over Th : %d, Total : %d\n",
		over_th_subgraphs_num,new_conc_list->size);

	*avg_diff_p=avg_diff;
	*max_diff_p=max_diff;
}



/********************************************************
* function : count_subgraphs_by_n_samples
*   search motifs probabilistic approach
*   use 'samples_num' as num of samples
*	results are returned in res list
* arguments:
*	N: Network to be searched
*	mtf_sz: motif size
*   samples_num :no of samples to sample
*	res : list pointer of results
*	net_type : REAL_NET - if real network
*              RAND_NET - if random network
* return values:
*	RC_OK - if not error occured
*   RC_ERR - if error occured
*********************************************************/
int
count_subgraphs_by_n_samples(Network *N, int mtf_sz, int samples_num, list64 **res_p, int net_type)
{
	int rc=RC_OK;
	int s,t,i,k;
	list *vrtx_set;
	list64 *res=*res_p;
	list_init(&vrtx_set);


	//find random connected subgraph of size mtf_sz
	//do it num of prob_total_ops
	for(i=1;i<=samples_num;i++) {
		//get a random edge to start
			if(DEBUG_LEVEL >=5 && GNRL_ST.out_log_flag)
				fprintf(GNRL_ST.log_fp,"sample #%d:\n",i);

			k=get_rand(N->edges_num);
			s=N->e_arr[k].s;
			t=N->e_arr[k].t;
			//self edge then pick random edge again
			if(s==t){
				i--;
				continue;
			}
			if(DEBUG_LEVEL >=6 && GNRL_ST.out_log_flag)
				fprintf(GNRL_ST.log_fp,"edge #1: (%d,%d)\n",s,t);
			//list_init(&ngbr_edges);
			//insert the two vertices of the edge to vertices set and call search subset
			list_insert(vrtx_set, s, NULL);
			list_insert(vrtx_set, t, NULL);
			if(GNRL_ST.efc_prob)
				rc = search_subset_prob_efc(N, vrtx_set, mtf_sz, res, net_type);
			else
				rc = search_subset_prob(N, vrtx_set, mtf_sz, res, net_type);
			//free hash history table
			//now remove these vertices from vertices set
			list_delete(vrtx_set, s);
			list_delete(vrtx_set, t);
			//list_free_mem(ngbr_edges);
			//if no big enough subset found then try again
			if (rc==RC_ERR)
				i--;
		}
	list_free_mem(vrtx_set);

	return RC_OK;
}


/********************************************************
* function : count_subgraphs_prob
*   search motifs probabilistic approach
*	results are returned in res list
* arguments:
*	N: Network to be searched
*	mtf_sz: motif size
*	res_p : refernce to list pointer of results
*	net_type : REAL_NET - if real network
*              RAND_NET - if random network
* return values:
*	RC_OK - if not error occured
*   RC_ERR - if error occured
*********************************************************/
int
count_subgraphs_prob(Network *N, int mtf_sz, list64 **res_p, int net_type, int rand_net_indx)
{
  //int rc=RC_OK;
	int curr_iter_samples_num; //current iteration samples num (convergness mode)
	int total_samples_num=0;  //total amples num (convergness mode)
	list64 *old_conc_list,*new_conc_list;
	double avg_diff,max_diff;
	list64 *res_list;


	//if not convergness mode
	if(GNRL_ST.prob_converge_mode==FALSE){
		count_subgraphs_by_n_samples(N,mtf_sz,GNRL_ST.prob_base_samples_num,res_p,net_type);
		GNRL_ST.prob_total_samples_num[rand_net_indx]=GNRL_ST.prob_base_samples_num;
	}else{
		//convergness mode

		//init old conc lists, in the beginnig its empty
		//this is equvalent to conc vec with all entries zero
		list64_init(&old_conc_list);
		curr_iter_samples_num=GNRL_ST.prob_base_samples_num;
		while(1) {
			list64_init(&res_list);
			count_subgraphs_by_n_samples(N,mtf_sz,curr_iter_samples_num,&res_list,net_type);

			total_samples_num+=curr_iter_samples_num;
			//NADAV : SHould be
			//GNRL_ST.prob_total_samples_num[rand_net_indx]=toal_samples_num;
			GNRL_ST.prob_total_samples_num[rand_net_indx]+=curr_iter_samples_num;
			//need to call join res now to get the concentration of each subgraph
			//it does not matter that I call it again in upper level cause then
			//it will have no work (need to check this works)
			join_subgraphs_res(&res_list,mtf_sz,rand_net_indx);
			//update the global res_tbl,if not first iteration then it already
			//hold the results till now (with all the previous iteration samples)
			update_global_res_tbl(*res_p,res_list,rand_net_indx);
			list64_free_mem(res_list);
			//update new conc list
			build_new_conc_list(&new_conc_list,*res_p);
			//calc diff between new and old conc lists
			calc_diff_conc_lists(new_conc_list,old_conc_list,&avg_diff,&max_diff,*res_p);

			//Two conditions to stop the iterations
			//
			//1) the average diff between the two vectors (lists) <=CONVERGNESS_DIFF_CONST
			//2) the max diff <= 2 * CONVERGNESS_DIFF_CONST
			//then exit the while
			if ((avg_diff<=GNRL_ST.prob_conv_diff) &&
				(max_diff<=(double)2*GNRL_ST.prob_conv_diff)){
				list64_free_mem(old_conc_list);
				list64_free_mem(new_conc_list);
				break;
			}else{
				//else raise num samples for the next iteration, and old conc list of next
				// iteration is the new of this iteration
				curr_iter_samples_num*=2;
				list64_free_mem(old_conc_list);
				old_conc_list=new_conc_list;
			}
			fprintf(stdout,"Avg diff = %.4f Max diff = %.4f ,Next sampling : %d\n",avg_diff,max_diff,curr_iter_samples_num);
		}
		fprintf(stdout,"Avg diff = %.4f Max diff = %.4f ,Total samples was %d\n",avg_diff,max_diff,total_samples_num);
	}
	return RC_OK;
}
