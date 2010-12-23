/************************************************************************
*
*  File name: main.c
*
*  Description: main file
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
#include "results.h"
#include "output.h"
#include "role.h"
#include "prob.h"
#include "random.h"
#include "switches.h"
#include "stubs.h"
#include "grassberger.h"
#include "clustering.h"


/*************************** Global variables ****************************/

int DEBUG_LEVEL=-50;


Network *G_N;
Network *G_WN;
Gnrl_st GNRL_ST;
//result table
Res_tbl RES_TBL;

char *input_network_fname;

//Gneralized motifs res table
Res_tbl GMTF_RES_TBL,GMTF_CMPLX_RES_TBL;

FILE *out_fp,*log_fp;
time_t start_time, end_time;

list64 *final_res, *final_res_all;

list64 *res_sub_motif;


//for grassberger
int died_colonies=0;
int over_populated_colonies=0;
double *weights_arr;



/******************************* Externs *********************************/
extern Hash *Role_hash;
extern role_st *Role_trans;
//Roles result table
extern Res_vec_tbl Roles_res_tbl;
//final results for role statistics
extern Role_res *Roles_final_res;
//Roles members 2 dim array [Role_id][node]
list **Roles_members;

//Gneralized motifs final res lists
extern list64 *gmtf_final_res;
extern list64 *gmtf_cmplx_final_res;

extern void
init_random_seed();

// return random integer in the interval 1 to max_val
extern int
get_rand(int max_val);
extern double get_rand_double();




/******************************* Declerations ****************************/
void
init_res_tbl(Res_tbl *res);
int
count_subgraphs_size_n(Network *N, int mtf_sz,list64 **res_p,int net_type,int rand_net_indx);
int
calc_final_results(Res_tbl*res_tbl, list64 **final_res_p, list64 **final_res_all_p, int rnd_net_num);
int
update_network(Network *RN,Network *N);


/******************************* Functions *******************************/


//start time meajuring.
//get as an argument a pointer to time_t where to store the start time
void
time_measure_start(Runtime *runtime)
{

	time(&runtime->start);
}

void
time_measure_stop(Runtime *runtime)
{
	time(&runtime->end);
	runtime->elapsed = difftime( runtime->end, runtime->start);
}


// free Network structure related allocated memory
void
free_subset_matrix_mem(Network *N)
{
	if(N!=NULL) {
		if(N->mat !=NULL){
			MatFree(N->mat);
			free(N->mat);
		}
	}
}

// free Network structure related allocated memory
void
free_network_mem(Network *N)
{
	int i;

	if(N!=NULL) {
		if(N->mat !=NULL){
			MatFree(N->mat);
			free(N->mat);
		}
		if(N->e_arr!=NULL)
			free(N->e_arr);
		if(N->e_arr_sin!=NULL)
			free(N->e_arr_sin);
		if(N->e_arr_dbl!=NULL)
			free(N->e_arr_dbl);
		if(N->indeg!=NULL)
			free(N->indeg);
		if(N->outdeg!=NULL)
			free(N->outdeg);
		if(N->doubledeg!=NULL)
			free(N->doubledeg);
		if(N->cluster_degree!=NULL)
			free(N->cluster_degree);
		if(N->roles_vec!=NULL) {
			for(i=1;i<N->vertices_num;i++)
				free(N->roles_vec[i]);
			free(N->roles_vec);
		}
		if(GNRL_ST.efc_prob) {
			if(N->e_map!=NULL){
				MatFree(N->e_map);
				free(N->e_map);
			}
			if(N->hub_edges_indices!=NULL)
				free(N->hub_edges_indices);
			if(N->hub_edges_vec!=NULL)
				free(N->hub_edges_vec);
		}
		if(GNRL_ST.r_conserve_layers){
			free(N->layer_edges_num);
			free(N->layer_node_num);
		}
	}
}


void
list_free_members_list(list *l)
{
	list_item *item;
	int *members_arr;
	for(item=list_get_next(l,NULL);item!=NULL;item=list_get_next(l,item)) {
		if(item->p !=NULL) {
			members_arr=(int*)item->p;
			free(members_arr);
			item->p=NULL;
		}
	}
}

void
list_free_all_members_list(list *l)
{
	list_item *item;
	Member *members_arr;
	for(item=list_get_next(l,NULL);item!=NULL;item=list_get_next(l,item)) {
		if(item->p !=NULL) {
			members_arr=(Member*)item->p;
			free(members_arr);
			item->p=NULL;
		}
	}
}

void
free_motif_st(Motif *mtf)
{
	if(mtf->members !=NULL) {
		list_free_mem(mtf->members);
		free(mtf->members);
	}
	if(mtf->all_members !=NULL) {
		list_free_mem(mtf->all_members);
		free(mtf->members);
	}
}








int*
create_subset_arr(list *org)
{
	int i;
	int *arr;
	list_item *org_item=NULL;

	if (org != NULL) {
		arr=(int*)calloc(org->size,sizeof(int));

		for(i=0,org_item=list_get_next(org,NULL); i<org->size;i++,org_item = list_get_next(org,org_item)) {
			arr[i]=org_item->val;
			if(arr[i]<0)
				at_exit(-1);
		}
		return arr;
	} else {
		return NULL;
	}
}



int64
get_sn_motif_id(Network *SN)
{
	int i;
	int mat_cells_num=SN->vertices_num*SN->vertices_num;
	int64 id=0;

	//id is the binary of the sub-matrix seen as one long array
	for(i=0;i<mat_cells_num;i++)
		id+=(int64)pow(2,i)*(int64)(SN->mat->full->m[i]);
	return id;
}


//check if equivalent set is in the hash
// important assumption - lists are sorted in increasing order
//returns : TRUE - if yes
//			FALSE - if not
int
check_hist_hash(Hash *hash, int *set_arr)
{
	int key,i;
	list_item *tmp;
	int *hist_set_arr;
#ifdef ENABLE_PROFILE
	calls_to_check_hist_hash++;
#endif
	if (set_arr==NULL)
		return FALSE;

	key=hash_get_key(hash,set_arr);
	if (hash->h[key] == NULL)
		return FALSE;

	for( tmp = list_get_next(hash->h[key],NULL); tmp!=NULL;
		tmp = list_get_next(hash->h[key],tmp) ) {
			hist_set_arr=(int*)tmp->p;
			//compare the lists
			//if all item vals are equal then return TRUE
			//otherwise return FALSE
			for(i=0;i<hash->bckt_sz;  i++) {
					if(set_arr[i] != hist_set_arr[i])
						break;
				}
			//if reached end of arrays - means they are identical
			if( i==hash->bckt_sz )
				return TRUE;
		}
	return FALSE;
}


//check if unique ,path through all  memebers of motif
//returns : TRUE - if unique
//			FALSE - otherwise
int
check_if_unique(Motif *mtf,int *mtf_vrtx_arr,int mtf_sz)
{
	int i,j;
	list_item *item;
	int*tmp_arr;
	int unique=TRUE;

	for(item=list_get_next(mtf->members,NULL);item!=NULL;
		item=list_get_next(mtf->members,item)) {
			tmp_arr=(int*)item->p;
			if (DEBUG_LEVEL >= 1 && GNRL_ST.out_log_flag) {
				fprintf(GNRL_ST.log_fp,"\n");
				for (i=0;i<mtf_sz;i++)
					fprintf(GNRL_ST.log_fp,"%d ",tmp_arr[i]);
				fprintf(GNRL_ST.log_fp,"\n");
			}
			for(i=0;i<mtf_sz && unique==TRUE;i++) {
				for(j=0;j<mtf_sz && unique==TRUE;j++) {
					if (tmp_arr[j]==mtf_vrtx_arr[i]) {
						unique=FALSE;
						break;
					}
				}
			}
		}
	return unique;
}


int
check_if_new_members(Motif *mtf,Member *mtf_vrtx_arr,int mtf_sz)
{
	int i;
	list_item *item;
	Member*tmp_arr;
	int identical=0;

	for(item=list_get_next(mtf->all_members,NULL);item!=NULL;
		item=list_get_next(mtf->all_members,item)) {
			tmp_arr=(Member*)item->p;
			if (DEBUG_LEVEL >= 1 && GNRL_ST.out_log_flag) {
				fprintf(GNRL_ST.log_fp,"\n");
				for (i=0;i<mtf_sz;i++)
					fprintf(GNRL_ST.log_fp,"%d ",tmp_arr[i].node);
				fprintf(GNRL_ST.log_fp,"\n");
			}
			identical=0;
			for(i=0;i<mtf_sz;i++) {
				if (tmp_arr[i].node==mtf_vrtx_arr[i].node) {
					identical++;
				}
			}
			//found identical members
			if (identical==mtf_sz)
				return FALSE;
		}
		//didnt find identical members
	return TRUE;
}


int64
get_rep_mtf_id(int64 mtf_id, int mtf_sz)
{
	int64 rep_mtf_id;
	list64 *iso_list;
	//get rep id
	switch(mtf_sz){
	case 3:
		rep_mtf_id=get_rep_mtf_id_3(mtf_id);
		break;
	case 4:
	case 5:
	case 6:
	case 7:
	case 8:
		iso_list=calc_mtf_id_iso(mtf_id, mtf_sz);
		rep_mtf_id=iso_list->l->val;
		break;
	default:
		rep_mtf_id=0;
		break;
	}
	return rep_mtf_id;
}


/********************************************************
* function : search_subset
*   Recursive function, search for subsets connected to the current vertex set
*   if vertex set equals motif size then end of recursion
*	results are returned in res list
* arguments:
*	N: Network to be searched
*	vrtx_set: list of vertices already participate in this search tree
*   hist_hash : hash table includes all sets of vertices already searched
*   mtf_sz: motif size
*	res : list pointer of results
*	net_type : REAL_NET - if real network
*              RAND_NET - if random network
* return values:
*	RC_OK - if not error occured
*   RC_ERR - if error occured
*********************************************************/
int
search_subset(Network *N, list *vrtx_set, Hash *hist_hash, int mtf_sz, list64 *res, int net_type)
{
	int i,k;
	list_item *item,*l_edge;
	list64_item *l_res_id;
	Network *SN;
	Motif *mtf;
	int64 mtf_id,rep_mtf_id=0;
	int role_id;
	Hash *subset_hash;
	int *cur_subset_arr;
	int *mtf_vrtx_arr;  //array that store the motif vertices
	Member *mtf_members_arr;
	int unique;
	int new_members;
	role_st *role=NULL;

	//number of vertices in set is equal to motif size
	//this is the recursive base
	if (mtf_sz == vrtx_set->size) {

		gen_subset_matrix(N, &SN, vrtx_set, mtf_sz);
		//calculate motif id
		mtf_id = get_sn_motif_id(SN);
		//update res tbl
		if ( (l_res_id=list64_get(res,(int64)mtf_id)) == NULL ){
			mtf=(Motif*)calloc(1,sizeof(Motif));
			mtf->id = mtf_id;
			mtf->count=1;
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

				//if omem but not for a specific subgraph
				//or if for specific subgrpah and this is the subgraph
				if( (GNRL_ST.specific_subgraph_members == 0) || (get_rep_mtf_id(mtf_id,mtf_sz)==GNRL_ST.specific_subgraph_members) ) {


					//alloc mtf_members_arr copy to it the motif vertices from vrtx_set
					mtf_members_arr = (Member*)calloc(mtf_sz,sizeof(Member));
					//get role rules for that id
					if(mtf_sz<=3)
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
							}
						}
						else
							mtf_members_arr[i].role=1;

					}
					//limit members list size by default
					if((GNRL_ST.list_members == TRUE) && (net_type == REAL_NET) && (mtf->all_members->size <= GNRL_ST.max_members_list_sz) ) {
						list_insert(mtf->all_members,1,(void*)mtf_members_arr);
						if (DEBUG_LEVEL >= 1 && GNRL_ST.out_log_flag) {
							fprintf(GNRL_ST.log_fp, "search_subset: new id added to res list %lld\n", mtf_id);
							fprintf(GNRL_ST.log_fp,"inserting motif all members: ");
							for (i=0;i<mtf_sz;i++) {
								fprintf(GNRL_ST.log_fp,"%d ",mtf_members_arr[i].node);
								fprintf(GNRL_ST.log_fp,"(%d) ",mtf_members_arr[i].role);
							}
							fprintf(GNRL_ST.log_fp,"\n");
						}
					} else {
						free(mtf_members_arr);
					}
				}
			}
			//now after mtf is fixed insert it to res list (in the p field)
			list64_insert(res,(int64)mtf_id,(void*)mtf);
		}else {

			//increase count by one
			((Motif*)(l_res_id->p))->count+=1;
			//if calc unique option is on
			//just add to motif memebers list the new motif vertices group
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
					list_insert(mtf->members,(int)((Motif*)(l_res_id->p))->count,(void*)mtf_vrtx_arr);
				} else {
					free(mtf_vrtx_arr);
				}
			}
			if( ((GNRL_ST.calc_roles) || ((GNRL_ST.list_members)&&(net_type == REAL_NET)) ) && (mtf_sz==GNRL_ST.mtf_sz) ) {
				//if omem but not for a specific subgraph
				//or if for specific subgrpah and this is the subgraph
				if( (GNRL_ST.specific_subgraph_members == 0) || (get_rep_mtf_id(mtf_id,mtf_sz)==GNRL_ST.specific_subgraph_members) ) {
					//alloc mtf_members_arr copy to it the motif vertices from vrtx_set
					mtf_members_arr = (Member*)calloc(mtf_sz,sizeof(Member));

					//get role rules for that id
					if(mtf_sz<=3)
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

							}
						}
						else
							mtf_members_arr[i].role=1;
					}
					mtf=(Motif*)l_res_id->p;
					//limit members list size by default value
					if((GNRL_ST.list_members == TRUE) && (net_type == REAL_NET) && (mtf->all_members->size <= GNRL_ST.max_members_list_sz)) {

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
		}
		//free subset matrix
		free_subset_matrix_mem(SN);
		free((void*)SN);
		return RC_OK;
	}
	//there is  a unique hash table according to the subset size
	//(fits the recursive depth)
	subset_hash=&hist_hash[vrtx_set->size+1];

	item = vrtx_set->l;
	//find neighbors off all vertices in vrtx_set
	for(k=1; k<=vrtx_set->size; k++) {


		//sparse matrix
		//pass through row
		for(l_edge=list_get_next(N->mat->spr->m[item->val].to,NULL);l_edge!=NULL;
		l_edge=list_get_next(N->mat->spr->m[item->val].to,l_edge)) {
			//if a "to" neighbor of the current vertices set
			//and not in vrtx_set already
			if(list_get(vrtx_set,l_edge->val)==NULL) {
				list_insert(vrtx_set,l_edge->val,NULL);
				//check if a such a subset was already counted
				// if already counted - then skip
				// if not counted then add this neighbor to set and call search_subset
				cur_subset_arr=create_subset_arr(vrtx_set);
				if (check_hist_hash(subset_hash, cur_subset_arr)==FALSE) {
					hash_insert(subset_hash,cur_subset_arr,cur_subset_arr);
					//recursive call
					search_subset(N,vrtx_set,hist_hash,mtf_sz,res,net_type);
				} else {
					free(cur_subset_arr);
				}
				list_delete(vrtx_set,l_edge->val);
			}
		}
		//pass through column
		for(l_edge=list_get_next(N->mat->spr->m[item->val].from,NULL);l_edge!=NULL;
		l_edge=list_get_next(N->mat->spr->m[item->val].from,l_edge)) {
			//if a "from" neighbor of the current vertices set
			//and not in vrtx_set already
			if(list_get(vrtx_set,l_edge->val)==NULL) {
				list_insert(vrtx_set,l_edge->val,NULL);
				//check if a such a subset was already counted
				// if already counted - then skip
				// if not counted then add this neighbor to set and call search_subset
				cur_subset_arr=create_subset_arr(vrtx_set);
				if (check_hist_hash(subset_hash, cur_subset_arr)==FALSE) {
					hash_insert(subset_hash,cur_subset_arr,cur_subset_arr);
					//recursive call
					search_subset(N,vrtx_set,hist_hash,mtf_sz,res,net_type);
				} else {
					free(cur_subset_arr);
				}
				list_delete(vrtx_set,l_edge->val);
			}
		}


		item = item->next;
	}
	return RC_OK;
}

/********************************************************
* function : count_subgraphs
*	results are returned in res list
* arguments:
*	N: Network to be searched
*	mtf_sz: motif size
*	res : list pointer of results
*	net_type : REAL_NET - if real network
*              RAND_NET - if random network
* return values:
*	RC_OK - if not error occured
*   RC_ERR - if error occured
*********************************************************/
int
count_subgraphs(Network *N, int mtf_sz, list64 **res_p, int net_type)
{
	int s,t,cnt,i;
	list *vrtx_set;
	list64_item *l_id;
	list64 *res=*res_p;
	//history hash for subsets
	Hash *hist_hash;
	int64 dummy;
	list_init(&vrtx_set);


	//hist_hash=multi_hash_init(mtf_sz, N->vertices_num+1);
	//pass through all edges in N
	for(i=1;i<=N->edges_num;i++) {

			s=N->e_arr[i].s;
			t=N->e_arr[i].t;
			if(s != t) /* if s==t this is self edge and self edges are automatically
				in the subgraphs if they exist */
			{
				//insert the two vertices of the edge to vertices set and call search subset
				list_insert(vrtx_set, s, NULL);
				list_insert(vrtx_set, t, NULL);
				//init hash history table - holds all subsets that were already
				//searched , starting from this edge
				hist_hash=multi_hash_init(mtf_sz);

				//printf("going to search with : %d  %d", s, t);
				search_subset(N, vrtx_set, hist_hash, mtf_sz, res, net_type);
				//free hash history table
				multi_hash_free_mem(hist_hash, mtf_sz);
				//now remove these vertices from vertices set
				list_delete(vrtx_set, s);
				list_delete(vrtx_set, t);
			}
	}
	list_free_mem(vrtx_set);
	//multi_hash_free_mem(hist_hash, mtf_sz);
	//pass through res list
	//for each id to get the real occurances number in the network
	// divide its score by the number of edges it contains
	//(num of edges is actually number of ones in the binary ID number)
	for(l_id=list64_get_next(res,NULL); l_id!=NULL; l_id=list64_get_next(res,l_id)) {
		cnt=0;
		dummy=(int64)l_id->val;
		if (dummy!=0) {
			for(i=0;i<64;i++,dummy>>=1) {
				if (dummy & 1)
					cnt++;
			}
			cnt -= ((Motif *)(l_id->p))->numberOfSelfEdges;
			((Motif*)(l_id->p))->count=(double)(((Motif*)(l_id->p))->count)/cnt;
		}
	}

	return RC_OK;
}


/********************************************************
* function : load_network
*	Load network form input file
* arguments:
*   N_p - reference to network structure
*   network_fname - network file name
*
* return values:
*	RC_OK - if not error occured
*   RC_ERR - if error occured
*********************************************************/
int
load_network(Network **N_p, char* network_fname)
{
	int   rc=RC_OK,s,t,i,j,weight;
	list *e_list;
	list_item *l_e;
	int max_node=0;
	Edge *e;
	FILE *fp;
	Network *N;
	int e_idx;
	int self_edge_exists=FALSE;
	int self_edge_number=0;
	int matVal = 0;

	fp = fopen(network_fname,"rt");
	if(fp==NULL) {
		printf("\nError: Cannot open input file : %s\n\tF input file name and path\n",network_fname);
		at_exit(-1);
	}
	N=(Network*)calloc(1,sizeof(Network));
	N->name=network_fname;
	N->edges_num = 0;

	//build edge array
	list_init(&e_list);
	//SRC->TRG format
	if(GNRL_ST.input_net_format==SRC_TRG_FORMAT) {
		if(GNRL_ST.quiet_mode==FALSE)
			fprintf(stdout,"\tReading Network file in <Source,Target,Weight> Format\n");
		//read all edges one by one and insert each edge to edge list
		while ( (rc = fscanf( fp,"%d %d %d\n", &s, &t, &weight) != EOF) ) {
			//insert edge to edge list
			e=(Edge*)calloc(1,sizeof(Edge));
			e->s=s;
			e->t=t;
			if(GNRL_ST.calc_weights == TRUE)
			{
				e->weight=weight;
			}
			else
			{
				e->weight=1;
			}
			list_insert(e_list,0,(void*)e);
			if (s>max_node) max_node=s;
			if (t>max_node) max_node=t;
			N->edges_num++;
			if(s == t)
			{
				self_edge_number++;
			}
		}
	} else {
		//TRG->SRC fvormat
		if(GNRL_ST.quiet_mode==FALSE)
			fprintf(stdout,"\tReading Network file in <Target,Source,Weight> Format\n");
		//read all edges one by one and insert each edge to edge list
		while ( (rc = fscanf( fp,"%d %d %d\n", &t, &s, &weight) != EOF) ) {
			//insert edge to edge list
			e=(Edge*)calloc(1,sizeof(Edge));
			e->s=s;
			e->t=t;
			if(GNRL_ST.calc_weights == TRUE)
			{
				e->weight=weight;
			}
			else
			{
				e->weight=1;
			}
			list_insert(e_list,0,(void*)e);
			if (s>max_node) max_node=s;
			if (t>max_node) max_node=t;
			N->edges_num++;
			if(s == t)
			{
				self_edge_number++;
			}
		}
	}

	fclose(fp);
	//directed graph
	if(GNRL_ST.undirected_flag==FALSE) {
		//allocate e_arr fill it according to e_list and free e_list
		N->e_arr=(Edge*)calloc(N->edges_num+1,sizeof(Edge));

		for(i=N->edges_num,l_e=list_get_next(e_list,NULL);i>0;i--,l_e=list_get_next(e_list,l_e)){
			N->e_arr[i].s=((Edge*)(l_e->p))->s;
			N->e_arr[i].t=((Edge*)(l_e->p))->t;
			N->e_arr[i].weight=((Edge*)(l_e->p))->weight;
		}
	} else {
		N->edges_num = N->edges_num*2 - self_edge_number;
		N->e_arr=(Edge*)calloc((N->edges_num+1),sizeof(Edge));
		
		for(i=N->edges_num,l_e=list_get_next(e_list,NULL);i>0;i--,l_e=list_get_next(e_list,l_e))
		{
			N->e_arr[i].s=((Edge*)(l_e->p))->s;
			N->e_arr[i].t=((Edge*)(l_e->p))->t;
			N->e_arr[i].weight=((Edge*)(l_e->p))->weight;
			if(((Edge*)(l_e->p))->s != ((Edge*)(l_e->p))->t) /* not a self edge */
			{
				i--;
				N->e_arr[i].s=((Edge*)(l_e->p))->t;
				N->e_arr[i].t=((Edge*)(l_e->p))->s;
				N->e_arr[i].weight=((Edge*)(l_e->p))->weight;
			}
		}
	}
	list_free_mem(e_list);

	N->vertices_num=max_node;

	rc |= MatInit(&N->mat,N->vertices_num,SPARSE);
	if(rc == RC_ERR) {
		printf("Error : Memory allocation failure\n");
		at_exit(-1);
	}

	//assign matrix entries according to edges
	for(i=1; i<=N->edges_num; i++) {
		if(GNRL_ST.calc_weights == TRUE) {
			MatAsgn(N->mat,N->e_arr[i].s,N->e_arr[i].t,N->e_arr[i].weight);
		}else{
			if(MatGet(N->mat,N->e_arr[i].s,N->e_arr[i].t) != 0) {
				printf("Error : Duplicate appearance of the edge (%d,%d) \n\tfound in the Network\n ",
				N->e_arr[i].s,N->e_arr[i].t);
				at_exit(-1);
			}else{
				MatAsgn(N->mat,N->e_arr[i].s,N->e_arr[i].t,1);
			}
		}
	}
	//if efc_prob app then
	//create e_map matrix mapping of edges [s,t] to their index in e_arr
	if(GNRL_ST.efc_prob){
		rc |= MatInit(&N->e_map,N->vertices_num,SPARSE);
		if(rc == RC_ERR) {
			printf("Error : Memory allocation failure\n");
			at_exit(-1);
		}
		for(i=1; i<=N->edges_num; i++)

			MatAsgn(N->e_map,N->e_arr[i].s,N->e_arr[i].t,i);
	}
	//check there are no Self edges- if there are any then output them to screen and  stop
	if(GNRL_ST.calc_self_edges == FALSE){
		for(i=1; i<=N->vertices_num; i++) {
			if(MatGet(N->mat,i,i)==1) {
				fprintf(stdout,"Self edges exist in Input Network!!!\n");
				fprintf(stdout,"Self Edges : (%d,%d)\t",i,i);
				self_edge_exists=TRUE;
			}
		}
		if(self_edge_exists==TRUE)
			at_exit(-1);
	}

	//allocate and fill arrays of single edges and double edges
	N->e_arr_sin=(Edge*)calloc(N->edges_num+1,sizeof(Edge));
	N->e_arr_dbl=(Edge*)calloc((N->edges_num+1),sizeof(Edge));
	N->e_sin_num=0;
	N->e_dbl_num=0;
	N->roots_num=0;
	N->leafs_num=0;
	//allocate indeg and out deg arrays
	N->indeg=(int*)calloc(N->vertices_num+2,sizeof(int));
	N->outdeg=(int*)calloc(N->vertices_num+2,sizeof(int));
	N->doubledeg=(int*)calloc(N->vertices_num+2,sizeof(int));
	if(GNRL_ST.calc_self_edges == TRUE){
		N->self_edge=(int*)calloc(N->vertices_num+2,sizeof(int));
	}
	//actually matrix is sparse anyway now
	if(N->mat->type == SPARSE) {
		for(i=1;i<=N->vertices_num;i++) {
			for(j=1;j<=N->vertices_num;j++) {
				if( (matVal = MatGet(N->mat,i,j))) {
					//if an edge and is not self edge
					if(i != j){
						//if the twin edge exists
						if(MatGet(N->mat,j,i)){
							//not inserted yet
							if(j>i) {
								//double edge- this way always the twin pair has indexes 2x-1,2x
								N->e_arr_dbl[++N->e_dbl_num].s=i;
								N->e_arr_dbl[N->e_dbl_num].t=j;
								N->e_arr_dbl[N->e_dbl_num].weight=matVal;
								N->e_arr_dbl[++N->e_dbl_num].s=j;
								N->e_arr_dbl[N->e_dbl_num].t=i;
								N->e_arr_dbl[N->e_dbl_num].weight=matVal;
							}
						}
						else {
							//single edge
							N->e_arr_sin[++N->e_sin_num].s=i;
							N->e_arr_sin[N->e_sin_num].t=j;
							N->e_arr_sin[N->e_sin_num].weight=matVal;
						}
					}
					else //self edge
					{
						N->e_arr_sin[++N->e_sin_num].s=i;
						N->e_arr_sin[N->e_sin_num].t=j;
						N->e_arr_sin[N->e_sin_num].weight=matVal;
					}
				}
			}
		}

		//fill in deg and out deg arrays
		for(i=0;i<=N->vertices_num;i++) {
			N->indeg[i]=0;
			N->outdeg[i]=0;
			if(GNRL_ST.calc_self_edges == TRUE)
			{
				N->self_edge[i]=FALSE;
			}
		}
		for (i=1; i<=N->vertices_num;i++){
			if(N->mat->spr->m[i].to==NULL)
				N->outdeg[i]=0;
			else
				N->outdeg[i]=N->mat->spr->m[i].to->size;
			if(N->mat->spr->m[i].from==NULL)
				N->indeg[i]=0;
			else
				N->indeg[i]=N->mat->spr->m[i].from->size;
			if((N->mat->spr->m[i].self_edge == 1) && (GNRL_ST.calc_self_edges == TRUE))
			{
				N->self_edge[i] = TRUE;
			}
		}
	}

	//statistics and global info about the network
	N->con_vertices_num=0;
	N->hub_deg=0;
	N->hub=0;
	N->in_hub_deg=0;
	N->in_hub=0;
	N->out_hub_deg=0;
	N->out_hub=0;
	//calc total num of connected vertices
	//and find hub preferenced
	for(i=1; i<=N->vertices_num;i++){
		if( (N->indeg[i]!=0) || (N->outdeg[i]!=0) )
			N->con_vertices_num++;
		if( ((N->indeg[i] + N->outdeg[i]) > N->hub_deg) ){
			N->hub_deg=N->indeg[i] + N->outdeg[i];
			N->hub=i;
		}
		if( N->indeg[i] > N->in_hub_deg){
			N->in_hub_deg=N->indeg[i];
			N->in_hub=i;
		}
		if( N->outdeg[i] > N->out_hub_deg){
			N->out_hub_deg=N->outdeg[i];
			N->out_hub=i;
		}
	}

	//if calc roles then init array or roles vector for all nodes
	if(GNRL_ST.calc_roles==TRUE && GNRL_ST.mtf_sz==3) {
		if( (N->roles_vec=(char**)calloc(N->vertices_num+1,sizeof(char*))) ==NULL)
				return RC_ERR;
		for(i=1;i<=N->vertices_num;i++)
			if( (N->roles_vec[i]=(char*)calloc(TOTAL_ROLES_3_NUM+1,sizeof(char))) ==NULL)
				return RC_ERR;
	}
	//if use clusterring in random network genereation
	if(GNRL_ST.use_clustering){
		N->cluster_degree=(double*)calloc(N->vertices_num+2,sizeof(double));
		for(i=0;i<=N->vertices_num;i++) {
			N->cluster_degree[i]=0.0;
		}
		/* fill clustering series */
		clustering_series(N,NULL);
	}
	//if efc prob approach then init hub_edges_vec and hub_edges_indices
	if(GNRL_ST.efc_prob){
		N->hub_edges_vec=(int*)calloc(N->edges_num+1,sizeof(int));
		N->hub_edges_indices=(int*)calloc(N->hub_deg+1,sizeof(int));

		//fill hub vector and hub indices
		i=0;
		for(l_e=list_get_next(N->e_map->spr->m[N->hub].to,NULL);l_e!=NULL;
			l_e=list_get_next(N->e_map->spr->m[N->hub].to,l_e)) {
					e_idx=*(int*)l_e->p;
					N->hub_edges_indices[++i]=e_idx;
					N->hub_edges_vec[e_idx]=i;
		}
		for(l_e=list_get_next(N->e_map->spr->m[N->hub].from,NULL);l_e!=NULL;
			l_e=list_get_next(N->e_map->spr->m[N->hub].from,l_e)) {
					e_idx=*(int*)l_e->p;
					N->hub_edges_indices[++i]=e_idx;
					N->hub_edges_vec[e_idx]=i;
		}
		if(i!=N->hub_deg)
			printf("Error in Hub degree\n");
	}
	//if conserve layers in random netowrks
	if(GNRL_ST.r_conserve_layers==TRUE){
		N->num_of_layers=GNRL_ST.r_layers_num;
		N->layer_node_num=(int*)calloc(N->num_of_layers+1,sizeof(int));
		N->layer_edges_num=(int*)calloc(N->num_of_layers+1,sizeof(int));
		//copy num of nodes in each layer
		for(i=1;i<=N->num_of_layers;i++)
			N->layer_node_num[i]=GNRL_ST.r_layer_sz[i];
		//run through e_arr and fill N->layer_edges_num
		//this is required for the randomizing process - in order to swtich edges
		//between nodes in the same layers only
		max_node=0;
		j=1;
		for(i=1;i<=N->num_of_layers;i++){
			max_node+=N->layer_node_num[i];
			while( (N->e_arr[j].s<=max_node) && (j<=N->edges_num) ){
				N->layer_edges_num[i]++;
				j++;
			}
		}
		//sanity check - that num of edges in  layers sum to total num of edges
		j=0;
		for(i=1;i<=N->num_of_layers;i++)
			j+=N->layer_edges_num[i];
		if(j!=N->edges_num){
			printf("\nERROR: in '-rcl' flag, check layers info\n");
			at_exit(-1);
		}
	}
	if(DEBUG_LEVEL>11)
		dump_network(stdout,N);
	*N_p=N;
	return rc;
}



int
allocate_network(Network *SRC, Network **TRG_p, char *trg_name)
{
	int i,rc=RC_OK;
	Network *TRG;
	//list_item *l_e;
	TRG=(Network*)calloc(1,sizeof(Network));
	TRG->name=trg_name;
	TRG->vertices_num = SRC->vertices_num;
	TRG->edges_num = SRC->edges_num;
	TRG->e_sin_num=SRC->e_sin_num;
	TRG->e_dbl_num=SRC->e_dbl_num;
	MatInit(&TRG->mat, TRG->vertices_num,SPARSE);

	//copy edge arr from source network
	TRG->e_arr=(Edge*)calloc(TRG->edges_num+1,sizeof(Edge));

	TRG->e_arr_sin=(Edge*)calloc((unsigned int)TRG->e_sin_num+1,sizeof(Edge));
	TRG->e_arr_dbl=(Edge*)calloc((unsigned int)TRG->e_dbl_num+1,sizeof(Edge));
	TRG->indeg=(int*)calloc((unsigned int)TRG->vertices_num+2,sizeof(int));
	TRG->outdeg=(int*)calloc((unsigned int)TRG->vertices_num+2,sizeof(int));
	TRG->doubledeg=(int*)calloc((unsigned int)TRG->vertices_num+2,sizeof(int));

	for(i=1;i<=TRG->vertices_num+1;i++) {
		TRG->indeg[i]=SRC->indeg[i];
		TRG->outdeg[i]=SRC->outdeg[i];
		TRG->doubledeg[i]=SRC->doubledeg[i];
	}
	//if calc roles then init array or roles vector for all nodes
	if(GNRL_ST.calc_roles==TRUE && GNRL_ST.mtf_sz==3) {
		if( (TRG->roles_vec=(char**)calloc(TRG->vertices_num+1,sizeof(char*))) ==NULL)
				return RC_ERR;
		for(i=1;i<=TRG->vertices_num;i++)
			if( (TRG->roles_vec[i]=(char*)calloc(TOTAL_ROLES_3_NUM+1,sizeof(char))) ==NULL)
				return RC_ERR;
	}
	if(GNRL_ST.use_clustering){
		TRG->cluster_degree=(double*)calloc(TRG->vertices_num+2,sizeof(double));
		for(i=1;i<=TRG->vertices_num+1;i++) {
			TRG->cluster_degree[i]=SRC->cluster_degree[i];
		}
	}
	//dump_network(TRG);
	*TRG_p=TRG;
	return rc;
}



/********************************************************
* function : duplicate_network
*	Duplicate
* arguments:
*   SRC - network source
*   TRG - netowrk target
*   network_fname - target network name
*
* return values:
*	RC_OK - if not error occured
*   RC_ERR - if error occured
*********************************************************/
int
duplicate_network(Network *SRC, Network **TRG_p, char *trg_name)
{
	int i,j,rc=RC_OK;
	Network *TRG;
	list_item *l_e;
	TRG=(Network*)calloc(1,sizeof(Network));
	TRG->name=trg_name;
	TRG->vertices_num = SRC->vertices_num;
	TRG->edges_num = SRC->edges_num;
	MatInit(&TRG->mat, TRG->vertices_num,SPARSE);
	//FULL: TRG->m = (char*)malloc(TRG->vertices_num*TRG->vertices_num*sizeof(char));
	//init matrix according to input network
	if(SRC->mat->type==FULL) {
		for (i=1; i<=TRG->vertices_num; i++)
			for (j=1; j<=TRG->vertices_num; j++)
				MatAsgn(TRG->mat,i,j,MatGet(SRC->mat,i,j));
	} else {
		//sparse
		for (i=1; i<=TRG->vertices_num; i++) {
			for(l_e=list_get_next(SRC->mat->spr->m[i].to,NULL);l_e!=NULL;
				l_e=list_get_next(SRC->mat->spr->m[i].to,l_e))
					MatAsgn(TRG->mat,i,l_e->val,*(int*)l_e->p);
				TRG->mat->spr->m[i].self_edge = SRC->mat->spr->m[i].self_edge;
		}
	}

	//copy edge arr from source network
	TRG->e_arr=(Edge*)calloc(TRG->edges_num+1,sizeof(Edge));
	
	for(i=1;i<=TRG->edges_num;i++) {
		TRG->e_arr[i].s=SRC->e_arr[i].s;
		TRG->e_arr[i].t=SRC->e_arr[i].t;
		TRG->e_arr[i].weight=SRC->e_arr[i].weight;
	}
	//copy e_arr_sin and e_arr_dup
	TRG->e_sin_num=SRC->e_sin_num;
	TRG->e_dbl_num=SRC->e_dbl_num;
	TRG->e_arr_sin=(Edge*)calloc((unsigned int)TRG->e_sin_num+1,sizeof(Edge));
	TRG->e_arr_dbl=(Edge*)calloc((unsigned int)TRG->e_dbl_num+1,sizeof(Edge));

	for(i=1;i<=TRG->e_sin_num;i++) {
		TRG->e_arr_sin[i].s=SRC->e_arr_sin[i].s;
		TRG->e_arr_sin[i].t=SRC->e_arr_sin[i].t;
		TRG->e_arr_sin[i].weight=SRC->e_arr_sin[i].weight;
	}
	for(i=1;i<=TRG->e_dbl_num;i++) {
		TRG->e_arr_dbl[i].s=SRC->e_arr_dbl[i].s;
		TRG->e_arr_dbl[i].t=SRC->e_arr_dbl[i].t;
		TRG->e_arr_dbl[i].weight=SRC->e_arr_dbl[i].weight;

	}
	TRG->indeg=(int*)calloc((unsigned int)TRG->vertices_num+2,sizeof(int));
	TRG->outdeg=(int*)calloc((unsigned int)TRG->vertices_num+2,sizeof(int));
	TRG->doubledeg=(int*)calloc((unsigned int)TRG->vertices_num+2,sizeof(int));
	for(i=1;i<=TRG->vertices_num+1;i++) {
		TRG->indeg[i]=SRC->indeg[i];
		TRG->outdeg[i]=SRC->outdeg[i];
		TRG->doubledeg[i]=SRC->doubledeg[i];
	}
	if(GNRL_ST.use_clustering){
		TRG->cluster_degree=(double*)calloc(TRG->vertices_num+2,sizeof(double));
		for(i=1;i<=TRG->vertices_num+1;i++) {
			TRG->cluster_degree[i]=SRC->cluster_degree[i];
		}
	}
	//if calc roles then init array or roles vector for all nodes
	if(GNRL_ST.calc_roles==TRUE && GNRL_ST.mtf_sz==3) {
		if( (TRG->roles_vec=(char**)calloc(TRG->vertices_num+1,sizeof(char*))) ==NULL)
				return RC_ERR;
		for(i=1;i<=TRG->vertices_num;i++)
			if( (TRG->roles_vec[i]=(char*)calloc(TOTAL_ROLES_3_NUM+1,sizeof(char))) ==NULL)
				return RC_ERR;
	}

	//copy hub params
	TRG->con_vertices_num=SRC->con_vertices_num;
	TRG->hub_deg=SRC->hub_deg;
	TRG->hub=SRC->hub;
	if(GNRL_ST.efc_prob){
		rc |= MatInit(&TRG->e_map,TRG->vertices_num,SPARSE);
		for(i=1; i<=TRG->edges_num; i++)
			MatAsgn(TRG->e_map,TRG->e_arr[i].s,TRG->e_arr[i].t,i);

		//if efc prob approach then allocate hub_edges_vec and hub_edges_indices
		TRG->hub_edges_vec=(int*)calloc(TRG->edges_num+1,sizeof(int));
		TRG->hub_edges_indices=(int*)calloc(TRG->hub_deg+1,sizeof(int));
		memcpy(TRG->hub_edges_vec,SRC->hub_edges_vec,(TRG->edges_num+1)*sizeof(int));
		memcpy(TRG->hub_edges_indices,SRC->hub_edges_indices,(TRG->hub_deg+1)*sizeof(int));

	}
	//if conserve layers mode then copy layers info
	if(GNRL_ST.r_conserve_layers==TRUE){
		TRG->num_of_layers=SRC->num_of_layers;
		TRG->layers_range=SRC->layers_range;

		TRG->layer_node_num=(int*)calloc(TRG->num_of_layers+1,sizeof(int));
		TRG->layer_edges_num=(int*)calloc(TRG->num_of_layers+1,sizeof(int));
		for(i=1;i<=TRG->num_of_layers;i++){
			TRG->layer_node_num[i]=SRC->layer_node_num[i];
			TRG->layer_edges_num[i]=SRC->layer_edges_num[i];
		}
	}

	//dump_network(TRG);
	*TRG_p=TRG;
	return RC_OK;
}







int
update_network(Network *RN,Network *N)
{
	int rc=RC_OK;
	int i,j,k,e_idx;
	list_item *l_e;
	int tot_edges_num=0;

	if (GNRL_ST.efc_prob || GNRL_ST.r_grassberger || GNRL_ST.use_stubs_method) {
		//dump_network(RN);

		for(j=1;j<=RN->vertices_num;j++) {
			RN->indeg[j]=0;
			RN->outdeg[j]=0;
		}
		for (j=1; j<=RN->vertices_num;j++){
			for(k=1;k<=RN->vertices_num;k++){
				if(MatGet(RN->mat,j,k)>0){
					RN->outdeg[j]++;
					RN->indeg[k]++;
					tot_edges_num++;
				}
			}
		}
		if(tot_edges_num!=N->edges_num){
			fprintf(stdout,"Total number of edges was not conserved!\n");
			rc=RC_ERR;
		}
		//check all nodes degree in the real networks and the rand network
		//are equal
		for(j=1; j<=RN->vertices_num;j++){
			if(RN->indeg[j]!=N->indeg[j]) {
				fprintf(stdout,"single node statistics not conserved ! , in deg node %d\n",j);
				fprintf(stdout,"RN indeg %d N indeg %d\n",RN->indeg[j],N->indeg[j]);
				if(GNRL_ST.out_log_flag)
					fprintf(GNRL_ST.log_fp,"single node statistics not conserved ! , in deg node %d\n",j);
				rc=RC_ERR;
			}

			if(RN->outdeg[j]!=N->outdeg[j]) {
				fprintf(stdout,"single node statistics not conserved ! , out deg node %d\n",j);
				fprintf(stdout,"RN outdeg %d N outdeg %d\n",RN->outdeg[j],N->outdeg[j]);
				if(GNRL_ST.out_log_flag)
					fprintf(GNRL_ST.log_fp,"single node statistics not conserved ! , out deg node %d\n",j);
				rc=RC_ERR;
			}
		}
	}


	//if efc_prob then update hub params
	//and hub edge vector and indices
	RN->con_vertices_num=0;
	RN->hub_deg=0;
	RN->hub=0;
	//calc total num of connected vertices
	//and find hub preferenced
	for(i=1; i<=RN->vertices_num;i++){
		if( (RN->indeg[i]!=0) || (RN->outdeg[i]!=0) )
			RN->con_vertices_num++;
		if( ((RN->indeg[i] + RN->outdeg[i]) > RN->hub_deg) ){
			RN->hub_deg=RN->indeg[i] + RN->outdeg[i];
			RN->hub=i;
		}
	}

	//if efc prob approach then init hub_edges_vec and hub_edges_indices
	if(GNRL_ST.efc_prob){
		//free e_map matrix and then init it according to new edges
		MatFree(RN->e_map);
		rc |= MatInit(&RN->e_map,RN->vertices_num,SPARSE);
		for(i=1; i<=RN->edges_num; i++)
			MatAsgn(RN->e_map,RN->e_arr[i].s,RN->e_arr[i].t,i);

		//free vec and indices and allocate (faster way to set to zero)
		free(RN->hub_edges_vec);
		free(RN->hub_edges_indices);
		//if efc prob approach then allocate hub_edges_vec and hub_edges_indices
		RN->hub_edges_vec=(int*)calloc(RN->edges_num+1,sizeof(int));
		RN->hub_edges_indices=(int*)calloc(RN->hub_deg+1,sizeof(int));

		//fill hub vector and hub indices
		i=0;
		for(l_e=list_get_next(RN->e_map->spr->m[RN->hub].to,NULL);l_e!=NULL;
			l_e=list_get_next(RN->e_map->spr->m[RN->hub].to,l_e)) {
					e_idx=*(int*)l_e->p;
					RN->hub_edges_indices[++i]=e_idx;
					RN->hub_edges_vec[e_idx]=i;
		}
		for(l_e=list_get_next(RN->e_map->spr->m[RN->hub].from,NULL);l_e!=NULL;
			l_e=list_get_next(RN->e_map->spr->m[RN->hub].from,l_e)) {
					e_idx=*(int*)l_e->p;
					RN->hub_edges_indices[++i]=e_idx;
					RN->hub_edges_vec[e_idx]=i;
		}
		if(i!=RN->hub_deg)
			printf("Error in Hub degree\n");
	}
	return rc;
}


int
process_rand_networks(Res_tbl *res_tbl,int mtf_sz)
{
	int rc=RC_OK;
	int i,j,l;
	Network *RN;
	int time_display_flag=0;
	time_t strt_time, curr_time;
	double elapsed_time, estimated_time;
	Res_tbl met_res_tbl;
	int real_vec13[14];
	int success;
	int dispair=0;
	char fname[1024];
	char *index_str;
	double *switch_ratio_arr; //array to store successfull switches ratio num of each rand net

	//statistics for the switches method
	switch_ratio_arr=(double*)calloc(GNRL_ST.rnd_net_num+1,sizeof(double));

	// init random seed
	init_random_seed();
	//for run time of random networks processing
	time_measure_start(&GNRL_ST.rand_net_time);
	//for time left approximation to the user
	time(&strt_time);

	for (i=0; i<GNRL_ST.rnd_net_num; i++) {
		// same triadic census ensamble.
		if((GNRL_ST.use_metropolis==TRUE) && (GNRL_ST.mtf_sz>3))
		{
			for (j=1;j<=13;j++)
				real_vec13[j]=0;
			//metropolis check - compute the 13 element vec on real network
			if(GNRL_ST.out_metrop_mat_flag==TRUE) {
				fprintf(GNRL_ST.mat_metrop_fp,"%d",G_N->vertices_num);
				fprintf(GNRL_ST.mat_metrop_fp,"\nThe real network :\n");
				for (l=1;l<=G_N->edges_num;l++)
					fprintf(GNRL_ST.mat_metrop_fp,"%d %d %d\n",G_N->e_arr[l].t,G_N->e_arr[l].s,1);
				fprintf(GNRL_ST.mat_metrop_fp,"End of real network\n");
			}

			init_res_tbl(&met_res_tbl);
			met_motifs_search_real(G_N,&met_res_tbl,real_vec13);
			if(GNRL_ST.out_metrop_mat_flag==TRUE)
				printf("real ok\n");
			list64_free_mem(met_res_tbl.real);

			gen_rand_network_metrop(&RN,real_vec13);
		}
		else if(GNRL_ST.use_clustering==TRUE)
			gen_rand_network_metrop_clustering(&RN);

		else if (GNRL_ST.use_stubs_method==TRUE)
		{
			// Use the stubs method
			success=0;
			dispair=0;
			while ((success==FALSE)&&(dispair<GEN_NETWORK_DISPAIR))
			{
				// Use Stubs method to generate random networks
				rc=gen_rand_network_stubs_method(&RN);
				if(rc==RC_OK)
					success=TRUE;
				else
					success=FALSE;
				if (success==FALSE)
					dispair++;

			}
			if (dispair==GEN_NETWORK_DISPAIR) // give up try switch method
				printf("Reached dispair limit on stubs method\n");

		}else {

			//Else use the default : the switches method
			if ((mtf_sz==2) || (GNRL_ST.r_dont_conserve_mutuals==TRUE) || (GNRL_ST.r_conserve_layers==TRUE) )
				// Generate random network without conserving mutual edges statistics
				gen_rand_network_switches_method(&RN,&switch_ratio_arr[i]);
			else
				// Generate random network conserving mutual edges statistics
				gen_rand_network_switches_method_conserve_double_edges(&RN,&switch_ratio_arr[i]);
		}
		//update network params (rest of the work that was not done in gen_rand_network..)
		//update degree sequences from the mat and check
		//for consistency with the real net
		rc=update_network(RN,G_N);
		if(rc!=RC_OK){
			printf("Error: Random network didnt conserve degree distribution\n");
			at_exit(-1);
		}

		//search the random network for subgraphs size n
		count_subgraphs_size_n(RN, mtf_sz, &res_tbl->rand_arr[i], RAND_NET, i+1);
		join_subgraphs_res(&res_tbl->rand_arr[i], mtf_sz,i+1);
		//fprintf(GNRL_ST.out_fp,"\nSummary motif results for rand %d network\n",i);
		//dump_motifs_res(RES_TBL.rand_arr[i], "random network");

		//if calc roles then update roles_res_tbl before freeing the network RN
		//this is done by nodes role vector
		if(GNRL_ST.calc_roles && GNRL_ST.mtf_sz==3)
			update_roles_res(RN->roles_vec,RN->vertices_num,Roles_res_tbl.rand[i]);


		if(GNRL_ST.out_all_rand_networks){
			strcat(fname,"");
			strncpy(fname,GNRL_ST.out_fname, strlen(GNRL_ST.out_fname)-4);
			fname[strlen(GNRL_ST.out_fname)-8]='\0';
			strcat(fname,"_RAND_");
			index_str=my_itoa(i+1,10);
			strcat(fname,index_str);
			strcat(fname,".txt");

			output_network_to_text_file(RN, fname);
		}
		free_network_mem(RN);
		free(RN);

		//res_tbl_mem_free(&RES_TBL);
		time_measure_stop(&GNRL_ST.total_time);

		//if intermediate output are required
		if(GNRL_ST.out_intermediate==TRUE) {
			//calculate final results and dump them to the results file
			calc_final_results(&RES_TBL, &final_res, &final_res_all,i+1);
			time_measure_stop(&GNRL_ST.total_time);

			dump_final_res(GNRL_ST.inter_out_fp,final_res,final_res_all,i+1);
		}

		if(GNRL_ST.quiet_mode==FALSE)
			printf(".");
		time_display_flag++;
		//throw some hint about time left (do not through if submotif searching)
		if(time_display_flag==GNRL_ST.rnd_net_num/10 ) {
			time(&curr_time);
			elapsed_time = difftime( curr_time, strt_time );
			estimated_time = elapsed_time * (GNRL_ST.rnd_net_num-(i+1)) / (i+1);
			if (estimated_time > 18000) {
				//display message in hours
				estimated_time /=3600;
				if (GNRL_ST.quiet_mode==FALSE)
					printf("\n Estimated run time left : %6.0f hours.\n", estimated_time);
			}
			else if (estimated_time > 300) {
				//display message in minutes
				estimated_time /=60;
				if (GNRL_ST.quiet_mode==FALSE)
					printf("\n Estimated run time left : %6.0f minutes.\n", estimated_time);
			}
			else {
				//display message in seconds
				if (GNRL_ST.quiet_mode==FALSE)
					printf("\n Estimated run time left : %6.0f seconds.\n", estimated_time );
			}
			time_display_flag=0;
			printf("\n");
		}
	}
	time_measure_stop(&GNRL_ST.rand_net_time);
	printf("\n");

	//calc statistics
	calc_stat(switch_ratio_arr,GNRL_ST.rnd_net_num,&GNRL_ST.rnstat.switch_ratio_mean,&GNRL_ST.rnstat.switch_ratio_stddev);
	free(switch_ratio_arr);

	return rc;
}

int
process_rand_networks_grassberger(Res_tbl *res_tbl,int mtf_sz, double *weights_arr)
{
	int rc=RC_OK;
	int i,j;
	Network *RN;
	int time_display_flag=0;
	time_t strt_time, curr_time;
	double elapsed_time, estimated_time;
	int success;
	int dispair=0;
	char fname[1024];
	char *index_str;
	int rnd_net_cnt=0;
	double *clones_time_arr; //array to store the clones time of the current rand net
	double *clones_mean_time_arr;
	double *clones_num_arr;
	double weights_sum=0;
	double eff_rand_nets;
	double max_weight;
	double dummy;
	FILE *fpc;

	// init random seed
	init_random_seed();

	//for time left approximation to the user
	time(&strt_time);
	printf("Generating random networks using Grassberger algorithm\n");
	printf("Colony size : %d\n",GNRL_ST.r_grass_colony_sz);
	clones_num_arr=(double*)calloc(NUM_OF_COLONIES_4_STATS+1,sizeof(double));
	clones_time_arr=(double*)calloc(MAX_CLONES_NUM+1,sizeof(double));
	clones_mean_time_arr=(double*)calloc(MAX_CLONES_NUM+1,sizeof(double));

	//special for statistics
	//Temp special arrangement
	//for pritnting the clones arr
	fname[0]='\0';
	strncpy(fname,GNRL_ST.out_fname, strlen(GNRL_ST.out_fname)-4);
	fname[strlen(GNRL_ST.out_fname)-8]='\0';
	strcat(fname,"_WEIGHTS.txt");
	fpc=fopen(fname,"w");

	printf("Calculate statistics about cloning times\n");
	//generate NUM_OF_COLONIES_4_STATS colonies to collect statistics about cloning times
	for(i=0;i<NUM_OF_COLONIES_4_STATS;i++){
		success=FALSE;
		dispair=0;
		while(success==FALSE){
			rc=stats_rand_network_grassberger_method(&RN,clones_time_arr,&dummy);

			if(rc==RC_OK) {
				success=TRUE;
				//update the mean time arr
				for(j=1;j<=MAX_CLONES_NUM;j++){
					//update the time to be the minimum if not zero
					if ( !((clones_mean_time_arr[j]==0) && (clones_time_arr[j]==0)) ){
						if( (clones_mean_time_arr[j]==0) && (clones_time_arr[j]>0) )
							clones_mean_time_arr[j]=clones_time_arr[j];
						else if ( (clones_mean_time_arr[j]>0) && (clones_time_arr[j]>0) ){
							if( clones_mean_time_arr[j]>clones_time_arr[j] )
								clones_mean_time_arr[j]=clones_time_arr[j];
						}
					}
				}
				//clones_mean_time_arr[j]=(double)((clones_mean_time_arr[j]*i+clones_time_arr[j]) / (double)(i+1));
				free_network_mem(RN);
				free(RN);
			}
			else
				success=FALSE;
			if (success==FALSE)
			{
				//printf("dispair - discard\n");
				dispair++;
				//set clones_time_arr to zero
				for(j=1;j<=MAX_CLONES_NUM;j++)
					clones_time_arr[j]=0;

			}
			if (dispair==GEN_NETWORK_GRASSBERGER_DISPAIR){
				printf("Reached dispair limit on grassberger method\n");
				at_exit(-1);
			}
		}
	}
	free(clones_num_arr);
	clones_num_arr=(double*)calloc(GNRL_ST.rnd_net_num+1,sizeof(double));

	//for run time of random networks processing
	time_measure_start(&GNRL_ST.rand_net_time);
	//repeat genrating random networks
	//each iteration generate colony and count subgraphs of all its networks
	for(i=0;i<GNRL_ST.rnd_net_num;i++){
		printf("\nApplying Grassberger algorithm colony %d \n",i+1);
		success=FALSE;
		dispair=0;
		while(success==FALSE){
			rc=gen_rand_network_grassberger_method(&RN,clones_mean_time_arr,&clones_num_arr[i],&weights_arr[i]);
			if(rc==RC_OK)
				success=TRUE;
			else
				success=FALSE;
			if (success==FALSE)
			{
				//printf("dispair - discard\n");
				dispair++;
			}
			if (dispair==GEN_NETWORK_GRASSBERGER_DISPAIR){
				printf("Reached dispair limit on grassberger method\n");
				at_exit(-1);
			}
		}

		//update network params (rest of the work that was not done in gen_rand_network)
		//update degree sequences from the mat and check
		//for consistency with the real net
		//entry here starts from 1 therefore I use i+1
		rc=update_network(RN,G_N);
		if(rc!=RC_OK){
			printf("Error: Random network didnt conserve degree distribution\n");
			at_exit(-1);
		}

		//search for subgraph size n
		count_subgraphs_size_n(RN, mtf_sz, &res_tbl->rand_arr[i], RAND_NET, i+1);
		join_subgraphs_res(&res_tbl->rand_arr[i], mtf_sz,i+1);

		//if calc roles then update roles_res_tbl before freeing the network RN
		//this is done by nodes role vector
		if(GNRL_ST.calc_roles && GNRL_ST.mtf_sz==3)
			update_roles_res(RN->roles_vec,RN->vertices_num,Roles_res_tbl.rand[i]);

		if(GNRL_ST.out_all_rand_networks){
			//strcat(fname,"");
			for(j=0;j<1024;j++){
				fname[j]='\0';
				//index_str[j]='\0';
			}
			strncpy(fname,GNRL_ST.out_fname, strlen(GNRL_ST.out_fname)-4);
			fname[strlen(GNRL_ST.out_fname)-8]='\0';
			strcat(fname,"_RAND_");
			index_str=my_itoa(i+1,10);
			strcat(fname,index_str);
			strcat(fname,".txt");

			output_network_to_text_file(RN, fname);
		}

		free_network_mem(RN);
		free(RN);
		rnd_net_cnt++;

		//dump time estimation
		time_display_flag++;
		//throw some hint about time left (do not through if submotif searching)
		if(time_display_flag==GNRL_ST.rnd_net_num/10) {
			time(&curr_time);
			elapsed_time = difftime( curr_time, strt_time );
			estimated_time = elapsed_time * (GNRL_ST.rnd_net_num-(i+1)) / (i+1);
			if (estimated_time > 18000) {
				//display message in hours
				estimated_time /=3600;
				if (GNRL_ST.quiet_mode==FALSE)
					printf("\n Estimated run time left : %6.0f hours.\n", estimated_time);
			}
			else if (estimated_time > 300) {
				//display message in minutes
				estimated_time /=60;
				if (GNRL_ST.quiet_mode==FALSE)
					printf("\n Estimated run time left : %6.0f minutes.\n", estimated_time);
			}
			else {
				//display message in seconds
				if (GNRL_ST.quiet_mode==FALSE)
					printf("\n Estimated run time left : %6.0f seconds.\n", estimated_time );
			}
			time_display_flag=0;
			printf("\n");
		}
	}
	//normalize the weights
	for(i=0;i<GNRL_ST.rnd_net_num;i++)
		weights_sum+=weights_arr[i];
	for(i=0;i<GNRL_ST.rnd_net_num;i++)
		weights_arr[i]=weights_arr[i]/weights_sum;
	//check the normalized sum is 1
	weights_sum=0;
	for(i=0;i<GNRL_ST.rnd_net_num;i++)
		weights_sum+=weights_arr[i];

	//calc eff_rand num
	max_weight=0;
	eff_rand_nets=0;
	for(i=0;i<GNRL_ST.rnd_net_num;i++)
		if(weights_arr[i]>max_weight)
			max_weight=weights_arr[i];
	for(i=0;i<GNRL_ST.rnd_net_num;i++)
		eff_rand_nets+=(weights_arr[i]/max_weight);

	GNRL_ST.rnstat.grass_reff=eff_rand_nets;
	//fprintf(stdout,"\nEffective Random Networks number is :%.2f\n",eff_rand_nets);
	if(DEBUG_LEVEL==-50 && GNRL_ST.out_log_flag){
		fprintf(stdout,"Weights array :\n");
		fprintf(GNRL_ST.log_fp,"Weights array :\n");
		for(i=0;i<GNRL_ST.rnd_net_num;i++){
			fprintf(stdout,"%.8f\n",weights_arr[i]);
			fprintf(GNRL_ST.log_fp,"%.8f\n ",weights_arr[i]);
		}
		fprintf(stdout,"\n");
		fprintf(GNRL_ST.log_fp,"\n");
	}
	for(i=0;i<GNRL_ST.rnd_net_num;i++)
		fprintf(fpc,"%.10f\n",weights_arr[i]);

	fclose(fpc);

	fprintf(stdout,"\nClones timing stats:\n");
	for(i=1;i<=MAX_CLONES_NUM && clones_mean_time_arr[i]>0;i++){
		fprintf(stdout,"%d ",(int)clones_mean_time_arr[i]);
		//fprintf(GNRL_ST.log_fp,"%d ",clones_num_arr[i]);
	}
	//calc statistics
	calc_stat(clones_num_arr,GNRL_ST.rnd_net_num,&GNRL_ST.rnstat.clones_mean,&GNRL_ST.rnstat.clones_stddev);

	fprintf(stdout,"\n");
	free(clones_num_arr);
	free(clones_time_arr);
	free(clones_mean_time_arr);
	time_measure_stop(&GNRL_ST.rand_net_time);
	time_measure_stop(&GNRL_ST.total_time);
	printf("\n");

	return rc;
}





int
process_input_args(int argc, char *argv[])
{
	int rc=RC_OK;
	int j,i=2;
	char in_fname_no_suffix[80]="";

	if ( argc<2 ) {
		print_public_help();
		//print_private_help()

		return RC_ERR;
	} else if(!strcmp(argv[1],"-h") || !strcmp(argv[1],"-help")) {
			print_public_help();
			at_exit(0);
	} else if(!strcmp(argv[1],"-hh")) {
			print_private_help();
			at_exit(0);
	} else {


		GNRL_ST.rnd_net_num = 0;
		input_network_fname=argv[1];

		//init GNRL_ST accoridng to default then overide
		//according to command line arguments
		GNRL_ST.mtf_sz=DEFAULT_MTF_SZ;
		GNRL_ST.rnd_net_num=DEFAULT_RAND_NETWORK_NUM;

		GNRL_ST.t_init=T0;
		GNRL_ST.iteration_factor=ITERATION_F;
		GNRL_ST.e_thresh=ETHRESH;
		GNRL_ST.use_stubs_method=FALSE;

		GNRL_ST.long_out_flag=FALSE;
		GNRL_ST.calc_unique_flag=TRUE;
		GNRL_ST.calc_roles=FALSE;
		GNRL_ST.input_net_format=SRC_TRG_FORMAT;
		GNRL_ST.undirected_flag=FALSE;
		GNRL_ST.calc_self_edges=FALSE;
		GNRL_ST.calc_weights=FALSE;
		GNRL_ST.run_prob_app=FALSE;
		GNRL_ST.prob_base_samples_num=0;
		GNRL_ST.prob_converge_mode=FALSE;
		GNRL_ST.prob_conv_diff=CONVERGNESS_DIFF_CONST;
		GNRL_ST.prob_conv_conc_th=CONC_THRESHOLD;
		GNRL_ST.unique_th=UNIQUE_TH;
		GNRL_ST.force_unique_th=FALSE;
		GNRL_ST.mfactor_th=MFACTOR_TH;
		GNRL_ST.zfactor_th=ZSCORE_TH;
		GNRL_ST.pval_th=PVAL_TH;
		GNRL_ST.max_members_list_sz=1000;
		GNRL_ST.top_motifs_list_sz=TOP_MOTIF_LIST_SZ;
		GNRL_ST.out_log_flag=FALSE;
		GNRL_ST.out_s_mat_flag=FALSE;
		GNRL_ST.out_l_mat_flag=FALSE;
		GNRL_ST.out_s_c_mat_flag=FALSE;
		GNRL_ST.out_members=FALSE;
		GNRL_ST.out_rand_mat_flag=FALSE;
		GNRL_ST.out_roles=FALSE;
		GNRL_ST.out_clustering=FALSE;
		GNRL_ST.quiet_mode=FALSE;
		GNRL_ST.use_metropolis=FALSE;
		GNRL_ST.use_clustering=FALSE;
		GNRL_ST.dont_search_real=FALSE;
		GNRL_ST.out_intermediate=FALSE;
		GNRL_ST.out_non_dangling_motifs=FALSE;
		GNRL_ST.r_switch_factor=(double)R_SWITCH_FACTOR;
		GNRL_ST.r_global_switch_mode=FALSE;
		GNRL_ST.list_members=FALSE;
		GNRL_ST.specific_subgraph_members=0; //no specific subgraph members
		GNRL_ST.efc_prob=FALSE;
		GNRL_ST.out_all_rand_networks=FALSE;
		GNRL_ST.r_grassberger=FALSE;
		GNRL_ST.r_grass_colony_sz=0;
		GNRL_ST.r_grass_colony_max_population=MAX_COLONY_SZ_RATIO;
		GNRL_ST.r_dont_conserve_mutuals=FALSE;
		GNRL_ST.dont_die=FALSE;
		GNRL_ST.r_conserve_layers=FALSE;
		GNRL_ST.r_layers_num=0;
		//if(!strcmp(strrchr(input_network_fname,"."),".txt")){
		if(strlen(input_network_fname)<4) {
			printf("Error: Network input file should have .txt suffix\n");
			at_exit(-1);
		}

		strncpy(in_fname_no_suffix,input_network_fname, strlen(input_network_fname)-4);
		strcat(in_fname_no_suffix,"");
		strcpy(GNRL_ST.out_fname,in_fname_no_suffix);
		strcpy(GNRL_ST.log_fname,in_fname_no_suffix);
		strcpy(GNRL_ST.mat_s_fname,in_fname_no_suffix);
		strcpy(GNRL_ST.mat_l_fname,in_fname_no_suffix);
		strcpy(GNRL_ST.mat_metrop_fname,in_fname_no_suffix);
		strcpy(GNRL_ST.inter_out_fname,in_fname_no_suffix);
		strcpy(GNRL_ST.members_fname,in_fname_no_suffix);
		strcpy(GNRL_ST.roles_fname,in_fname_no_suffix);
		strcpy(GNRL_ST.rand_all_fname,in_fname_no_suffix);
		strcpy(GNRL_ST.clust_fname,in_fname_no_suffix);

		strcat(GNRL_ST.out_fname,"_OUT.txt");
		strcat(GNRL_ST.log_fname,"_LOG.txt");
		strcat(GNRL_ST.mat_s_fname,"_MAT.txt");
		strcat(GNRL_ST.mat_l_fname,"_MAT_long.txt");
		strcat(GNRL_ST.mat_metrop_fname,"_METROP.txt");
		strcat(GNRL_ST.inter_out_fname,"_INTER.txt");
		strcat(GNRL_ST.members_fname,"_MEMBERS.txt");
		strcat(GNRL_ST.roles_fname,"_ROLES.txt");
		strcat(GNRL_ST.rand_all_fname,"_MAT_R.txt");
		strcat(GNRL_ST.clust_fname,"_CLUST.txt");

		//if input file name is network_example then set uniqueness threshold to 4
		if(!strcmp(input_network_fname,"network_exmp.txt"))
			GNRL_ST.unique_th=4;

		//interpret command line arguements
		while(i<argc) {
			//motif size
			if (!strcmp(argv[i],"-s")) {
				if(++i>=argc) {
					printf("Error :need to correct -s <num> flag\n");
					return RC_ERR;
				}
				GNRL_ST.mtf_sz = atoi(argv[i]);
				if(GNRL_ST.mtf_sz > 8 || GNRL_ST.mtf_sz < 2) {
					printf("Error :This version of mfinder supports motif size 3 to 8 only\nChange argument after -s flag\n");
					return RC_ERR;
				}

				//random network num
			} else if (!strcmp(argv[i],"-r")) {
				if(++i>=argc) {
					printf("Error :need to correct -r <num> flag\n");
					return RC_ERR;
				}
				GNRL_ST.rnd_net_num = atoi(argv[i]);
				if(GNRL_ST.rnd_net_num < 0) {
					printf("Error :Number of random networks should be a positive number\nChange argument after -r flag\n");
					return RC_ERR;
				}
				// metropolis iteration numbers
			} else if (!strcmp(argv[i],"-iter")) {
				GNRL_ST.iteration_factor = atof(argv[++i]);
				// initial metropolis temperature
			} else if (!strcmp(argv[i],"-rs")) {
				GNRL_ST.use_stubs_method = TRUE;
			}else if (!strcmp(argv[i],"-t0")) {
				GNRL_ST.t_init = atof(argv[++i]);
				// metropolis energy threshold
			} else if (!strcmp(argv[i],"-eth")) {
				GNRL_ST.e_thresh = atof(argv[++i]);
				//unique th
			} else if (!strcmp(argv[i],"-u")) {
				if(GNRL_ST.calc_unique_flag==TRUE) {
					GNRL_ST.unique_th=atoi(argv[++i]);
					GNRL_ST.force_unique_th=TRUE;
				}
				//dont calc uniqueness
			} else if (!strcmp(argv[i],"-nu")) {
				GNRL_ST.unique_th=0;
				GNRL_ST.calc_unique_flag=FALSE;
				//work with self edges
			} else if (!strcmp(argv[i],"-se")) {
				GNRL_ST.calc_self_edges=TRUE;
				//input format is target source
			} else if (!strcmp(argv[i],"-ts")) {
				GNRL_ST.input_net_format=TRG_SRC_FORMAT;
				//work with weights
			} else if (!strcmp(argv[i],"-w")) {
				GNRL_ST.calc_weights=TRUE;
				//zscore th to use
			} else if (!strcmp(argv[i],"-z")) {
				GNRL_ST.zfactor_th = atof(argv[++i]);
				//mfactor th to use
			} else if (!strcmp(argv[i],"-m")) {
				GNRL_ST.mfactor_th = atof(argv[++i]);
				//undirected graph
			} else if (!strcmp(argv[i],"-nd")) {
				GNRL_ST.undirected_flag=TRUE;
			//probabilistic approach
			}else if (!strcmp(argv[i],"-cr")) {
				if(GNRL_ST.mtf_sz==3) {
					GNRL_ST.calc_roles=TRUE;
					GNRL_ST.out_roles=TRUE;
				}else{
					printf("Error: Cannot use -cr flag for more then 3 nodes motifs\n");
					return RC_ERR;
				}
			//probabilistic approach
			} else if (!strcmp(argv[i],"-pold")) {
				GNRL_ST.run_prob_app=TRUE;
				GNRL_ST.prob_base_samples_num = atoi(argv[++i]);
				//efficient samples method
			} else if (!strcmp(argv[i],"-p")) {
				GNRL_ST.run_prob_app=TRUE;
				GNRL_ST.efc_prob=TRUE;
				GNRL_ST.prob_base_samples_num = atoi(argv[++i]);
				//if convergness probabilistic mode
			} else if (!strcmp(argv[i],"-pconv")) {
				GNRL_ST.prob_converge_mode=TRUE;
				//converged conc th
			} else if (!strcmp(argv[i],"-pconv_th")) {
				GNRL_ST.prob_conv_conc_th = atof(argv[++i]);
				//convergeness diff requirement
			} else if (!strcmp(argv[i],"-pconv_diff")) {
				GNRL_ST.prob_conv_diff = atof(argv[++i]);
				//quiet mode
			} else if (!strcmp(argv[i],"-q")) {
				GNRL_ST.quiet_mode=TRUE;
			//dont die mode
			} else if (!strcmp(argv[i],"-dd")) {
				GNRL_ST.dont_die=TRUE;
			//out .MAT in conc not counts
			} else if (!strcmp(argv[i],"-oc")) {
				GNRL_ST.out_s_c_mat_flag=TRUE;
			} else if (!strcmp(argv[i],"-olog")) {
				GNRL_ST.out_log_flag=TRUE;
			} else if (!strcmp(argv[i],"-omat")) {
				GNRL_ST.out_s_mat_flag=TRUE;
			} else if (!strcmp(argv[i],"-orall")) {
				GNRL_ST.out_rand_mat_flag=TRUE;
			}else if (!strcmp(argv[i],"-omet")) {
				GNRL_ST.out_metrop_mat_flag=TRUE;
			} else if (!strcmp(argv[i],"-met")) {
				GNRL_ST.use_metropolis=TRUE;
			}else if (!strcmp(argv[i],"-rclust")) {
				GNRL_ST.use_clustering=TRUE;
			}else if (!strcmp(argv[i],"-oclust")) {
				GNRL_ST.out_clustering=TRUE;
				//output all random networks
			} else if (!strcmp(argv[i],"-ornet")) {
				GNRL_ST.out_all_rand_networks=TRUE;
				//output members file - list all members of all subgraphs
			} else if (!strcmp(argv[i],"-omem")) {
				GNRL_ST.list_members=TRUE;
				GNRL_ST.out_members=TRUE;
				//list members of a specific subgraph only
			} else if (!strcmp(argv[i],"-ospmem")) {
				GNRL_ST.list_members=TRUE;
				GNRL_ST.out_members=TRUE;
				GNRL_ST.specific_subgraph_members=atoi(argv[++i]);
				//dont search on real network flag
			} else if (!strcmp(argv[i],"-nor")) {
				GNRL_ST.dont_search_real=TRUE;
				//output intermediate results
			} else if (!strcmp(argv[i],"-oi")) {
				GNRL_ST.out_intermediate=TRUE;
				//output list of all deteced non-dangling motifs
			} else if (!strcmp(argv[i],"-onodangl")) {
				GNRL_ST.out_non_dangling_motifs=TRUE;
				//top motifs lists size
			} else if (!strcmp(argv[i],"-otop")) {
				GNRL_ST.top_motifs_list_sz = atoi(argv[++i]);
				//unlimitedMax members list size
			} else if (!strcmp(argv[i],"-maxmem")) {
				GNRL_ST.max_members_list_sz = atoi(argv[++i]);
				//Global switch mode inrandomizing.
			} else if (!strcmp(argv[i],"-nsr")) {
				GNRL_ST.r_switch_factor = atof(argv[++i]);
				//Global switch mode inrandomizing.
			} else if (!strcmp(argv[i],"-rgs")) {
				GNRL_ST.r_global_switch_mode = TRUE;
				//Grassberger alg for generting random networks
			} else if (!strcmp(argv[i],"-rgrass")) {
				GNRL_ST.r_grassberger=TRUE;
				GNRL_ST.r_grass_colony_sz = atoi(argv[++i]);
			} else if (!strcmp(argv[i],"-rgrass_max_sz")) {
				GNRL_ST.r_grass_colony_max_population=atoi(argv[++i]);
			} else if (!strcmp(argv[i],"-rdm")) {
				GNRL_ST.r_dont_conserve_mutuals=TRUE;
			//conserve layers in random networks
			} else if (!strcmp(argv[i],"-rcl")) {
				GNRL_ST.r_conserve_layers=TRUE;
				GNRL_ST.r_layers_num=atoi(argv[++i]);
				GNRL_ST.r_layer_sz=(int*)calloc(GNRL_ST.r_layers_num+1,sizeof(int));
				for(j=1;j<=GNRL_ST.r_layers_num;j++)
					GNRL_ST.r_layer_sz[j]=atoi(argv[++i]);
			}else if (!strcmp(argv[i],"-f")) {
				if(++i>=argc) {
					printf("Error :need to correct -f <output file name> flag\n");
					return RC_ERR;
				}
				strcpy(GNRL_ST.out_fname,argv[i]);
				strcpy(GNRL_ST.log_fname,argv[i]);
				strcpy(GNRL_ST.mat_s_fname,argv[i]);
				strcpy(GNRL_ST.mat_l_fname,argv[i]);
				strcpy(GNRL_ST.mat_metrop_fname,argv[i]);
				strcpy(GNRL_ST.inter_out_fname,argv[i]);
				strcpy(GNRL_ST.members_fname,argv[i]);
				strcpy(GNRL_ST.roles_fname,argv[i]);
				strcpy(GNRL_ST.rand_all_fname,argv[i]);

				strcat(GNRL_ST.out_fname,"_OUT.txt");
				strcat(GNRL_ST.log_fname,"_LOG.txt");
				strcat(GNRL_ST.mat_s_fname,"_MAT.txt");
				strcat(GNRL_ST.mat_l_fname,"_MAT_long.txt");
				strcat(GNRL_ST.mat_metrop_fname,"_METROP.txt");
				strcat(GNRL_ST.inter_out_fname,"_INTER.txt");
				strcat(GNRL_ST.members_fname,"_MEMBERS.txt");
				strcat(GNRL_ST.roles_fname,"_ROLES.txt");
				strcat(GNRL_ST.rand_all_fname,"_MAT_R.txt");
			} else if (!strcmp(argv[i],"-ol")) {
				GNRL_ST.out_l_mat_flag=TRUE;
			} else {
				print_public_help();
				//print_private_help();
				at_exit(-1);
			}
			i++;
		}
	}

	if(GNRL_ST.quiet_mode==FALSE)
			printf("Input Network file is %s\n", input_network_fname);
	return rc;
}

int
gnrl_init()
{
	int rc=RC_OK;
	//open output  file
	GNRL_ST.out_fp=fopen(GNRL_ST.out_fname, "w");
	if(GNRL_ST.out_fp==NULL) {
		printf("\nError: Cannot open outpuf file: %s\n\t Check output file name and path\n",GNRL_ST.out_fname);
		at_exit(-1);
	}

	if(GNRL_ST.run_prob_app && GNRL_ST.out_s_mat_flag)
		GNRL_ST.out_s_c_mat_flag=TRUE;

	if(GNRL_ST.out_log_flag==TRUE)
		GNRL_ST.log_fp=fopen(GNRL_ST.log_fname, "w");
	if(GNRL_ST.out_s_mat_flag==TRUE)
		GNRL_ST.mat_s_fp=fopen(GNRL_ST.mat_s_fname, "w");
	if(GNRL_ST.out_l_mat_flag==TRUE)
		GNRL_ST.mat_l_fp=fopen(GNRL_ST.mat_l_fname, "w");
	if(GNRL_ST.out_metrop_mat_flag==TRUE)
		GNRL_ST.mat_metrop_fp=fopen(GNRL_ST.mat_metrop_fname, "w");
	if(GNRL_ST.out_intermediate==TRUE)
		GNRL_ST.inter_out_fp=fopen(GNRL_ST.inter_out_fname, "w");
	if(GNRL_ST.out_members==TRUE)
		GNRL_ST.members_fp=fopen(GNRL_ST.members_fname, "w");
	if(GNRL_ST.out_roles==TRUE)
		GNRL_ST.roles_fp=fopen(GNRL_ST.roles_fname, "w");
	if(GNRL_ST.out_rand_mat_flag==TRUE)
		GNRL_ST.mat_rand_fp=fopen(GNRL_ST.rand_all_fname, "w");
	if(GNRL_ST.out_clustering==TRUE)
		GNRL_ST.clust_fp=fopen(GNRL_ST.clust_fname, "w");

	init_res_tbl(&RES_TBL);


	if(GNRL_ST.run_prob_app)
		GNRL_ST.prob_total_samples_num=(int*)calloc(GNRL_ST.rnd_net_num+1,sizeof(int));

	if(GNRL_ST.list_members==TRUE) {
		//currently support only size 3
		init_role_hash(GNRL_ST.mtf_sz);
		init_role_trans(GNRL_ST.mtf_sz);
	}
	if(GNRL_ST.calc_roles==TRUE) {
		//currently support only size 3
		init_role_hash(GNRL_ST.mtf_sz);
		init_role_trans(GNRL_ST.mtf_sz);
		init_roles_res_tbl(GNRL_ST.rnd_net_num);
		init_roles_members(GNRL_ST.mtf_sz);
	}
	//if probabilistic method and uniqueness threshold was not set by command option
	//then use th=0;   in order not to miss motifs when run with small number of samples
	if( (GNRL_ST.run_prob_app) && (!GNRL_ST.force_unique_th) ) {
		GNRL_ST.unique_th=0;
	}

	return rc;
}



/********************************************************
* function : count_subgrpah_size_n
*   Search network for subgraphs size mtf_sz
*	results are returned in res list
* arguments:
*	N: Network to be searched
*	mtf_sz: motif size
*	res : list pointer of results
*	net_type : REAL_NET - if real network
*              RAND_NET - if random network
*   rand_net_indx : index of the current random network (starts from 1)
* return values:
*	RC_OK - if not error occured
*   RC_ERR - if error occured
*********************************************************/
int
count_subgraphs_size_n(Network *N, int mtf_sz,list64 **res_p,int net_type,int rand_net_indx)
{
	int rc=RC_OK;

	if(GNRL_ST.run_prob_app==FALSE)
		//exhaustive search
		count_subgraphs(N, mtf_sz, res_p, net_type);
	else
		//probabilistic alg (sampling)
		count_subgraphs_prob(N, mtf_sz, res_p, net_type,rand_net_indx);
	return rc;
}


//process subgraph search on Real network
//summarize results for real
int
motifs_search_real(Network *N)
{
	int rc=RC_OK;

	time_measure_start(&GNRL_ST.real_net_time);
	//search motif size n
	count_subgraphs_size_n(N, GNRL_ST.mtf_sz,&RES_TBL.real,REAL_NET,0);
	//calc result after isomorphism of ids
	join_subgraphs_res(&RES_TBL.real, GNRL_ST.mtf_sz, 0);

	//if calc roles option
	if(GNRL_ST.calc_roles && GNRL_ST.mtf_sz==3) {
			update_roles_res(N->roles_vec,N->vertices_num,Roles_res_tbl.real);
			rc |= fill_roles_members(N->vertices_num,N->roles_vec,Roles_res_tbl.real,Roles_members);
			Roles_final_res=(Role_res*)calloc(TOTAL_ROLES_3_NUM+1,sizeof(Role_res));
	}

	free_network_mem(N);
	free(N);

	fprintf(GNRL_ST.out_fp,"\n   Summary motif results\n");
	fprintf(GNRL_ST.out_fp,"   =====================\n");

	time_measure_stop(&GNRL_ST.real_net_time);
	if(GNRL_ST.quiet_mode==FALSE)
		dump_time_measure(stdout, "Real network processing runtime was:", &GNRL_ST.real_net_time);

	return rc;
}

char* my_itoa(int val, int base)
{
                static char buf[32] = {0};
                int i = 30;
                for(; val && i ; --i, val /= base)
                        buf[i] = "0123456789abcdef"[val % base];
                return &buf[i+1];
}


int
at_exit(int rc) {
	if (GNRL_ST.dont_die) {
		fprintf(stdout,"\nPress 'ENTER' to terminate\n");
			fgetc(stdin);
	}
	exit(rc);
}

int
main(int argc, char *argv[])
{
	int rc=0;
	Network *N;

	time_measure_start(&GNRL_ST.total_time);

	//process input arguments and init global structure
	rc |=process_input_args(argc,argv);

	if(GNRL_ST.quiet_mode==FALSE)
		printf("mfinder Version %.2f\n\n",VERSION);


	if(rc==RC_ERR)
		at_exit(-1);

	//general initialization
	rc|=gnrl_init();
	if(rc==RC_ERR)
		at_exit(-1);

	// load network from input file
	if(GNRL_ST.quiet_mode==FALSE)
		printf("Loading Network\n");
	load_network(&G_N,input_network_fname);
	duplicate_network(G_N,&N,"real_network");

	init_random_seed();

	if(rc==RC_ERR)
		at_exit(-1);
	if(GNRL_ST.quiet_mode==FALSE)
		printf("Searching motifs size %d\nProcessing Real network...\n",GNRL_ST.mtf_sz);

	//search motifs size n in Real network
	if(GNRL_ST.dont_search_real!=TRUE)
		rc|=motifs_search_real(N);
	if(rc==RC_ERR)
		at_exit(-1);


	if(GNRL_ST.quiet_mode==FALSE)
		printf("Processing Random networks\n");
	if (GNRL_ST.rnd_net_num>0) {
		// create random networks with same single node statisticfs as the input network
		if(GNRL_ST.r_grassberger==FALSE){
			//use switches or stubs
			rc|=process_rand_networks(&RES_TBL, GNRL_ST.mtf_sz);
		} else {
			//use grassberger alg
			weights_arr=(double*)calloc(GNRL_ST.rnd_net_num+1,sizeof(double));
			rc|=process_rand_networks_grassberger(&RES_TBL, GNRL_ST.mtf_sz,weights_arr);
		}
		if(rc==RC_ERR)
			at_exit(-1);
	}

	if(GNRL_ST.quiet_mode==FALSE)
		printf("Calculating Results...\n");
	if(GNRL_ST.rnd_net_num>=0){
		//calculate final results and dump them to the results file
		if(!GNRL_ST.r_grassberger) {
			calc_final_results(&RES_TBL, &final_res, &final_res_all,GNRL_ST.rnd_net_num);
		} else {
			//Nadav change for GRASS NEW
			calc_final_results_grassberger(&RES_TBL, FALSE, res_sub_motif, &final_res, &final_res_all,GNRL_ST.rnd_net_num,weights_arr);
		}
	}

	//calculate final results
	time_measure_stop(&GNRL_ST.total_time);
	//output results
	rc|=output_results(final_res,final_res_all);

	free_network_mem(G_N);
	free(G_N);

	final_res_free(final_res);
	final_res_free(final_res_all);

	if(GNRL_ST.r_grassberger)
		free(weights_arr);

	res_tbl_mem_free(&RES_TBL);

	if(GNRL_ST.calc_roles==TRUE) {
		free_mem_role_hash();
		free_roles_res_tbl(GNRL_ST.rnd_net_num);
		free_role_members();
	}

	if(rc==RC_ERR)
		at_exit(-1);

	at_exit(rc);
	//return rc;

}
