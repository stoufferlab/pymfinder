/************************************************************************
*
*  File name: metropolis.c
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


void zero_neighbours(int *arr,int *arr_len)
{
	register int i;

	for (i=0;i<ARR_NUM;i++)
	{
		arr[i]=0;
	}

	*arr_len=0;
}

void fill_neighbours(Network *N,int s,int t,int direct,int *arr,int *arr_len)
/* add all weakly connected nodes of s to arr, and return also arr_len */
/* direct - 1 if  there is an edge from s to t.
			0 if there is an edge from t to s
don't add node t */
{
	list_item *tmp;
	list_item *tmp2;

	// add all incoming and outgoing edges from/to node s, except t
	tmp = list_get_next(N->mat->spr->m[s].to,NULL);
	tmp2 = list_get_next(N->mat->spr->m[s].from,NULL);

	while (tmp != NULL)
	{
		if (direct)
		{
			if (tmp->val!=t)
			{
				arr[*arr_len]=tmp->val;
				(*arr_len)=(*arr_len)+1;
			}
		}
		else
		{
			arr[*arr_len]=tmp->val;
			(*arr_len)=(*arr_len)+1;

		}
		tmp=list_get_next(N->mat->spr->m[s].to,tmp);
	}

	while (tmp2 != NULL)
	{
		if (!direct)
		{
			if (tmp2->val!=t)
			{
				arr[*arr_len]=tmp2->val;
				(*arr_len)=(*arr_len)+1;
			}
		}
		else
		{
			arr[*arr_len]=tmp2->val;
			(*arr_len)=(*arr_len)+1;
		}
		tmp2=list_get_next(N->mat->spr->m[s].from,tmp2);
	}

}

//   computes the current 13 3-node circuits of s1-t1 and neighbours
//   and s2-t2+neighbours. Puts results in sub-a 13*1 vector.
//   computes the future 13 3-node circuits of s1-t2 and neighbours
//   and s2-t1+neighbours. Puts results in add-a 13*1 vector.
//   Takes into account that single edges are switched
void fill_13(Network *N,int s1,int t1,int s2,int t2,int sub[14],int add[14])

{
   int sarr1[ARR_NUM];
   int tarr1[ARR_NUM];
   int sarr2[ARR_NUM];
   int tarr2[ARR_NUM];
   int sarr1_len,tarr1_len,sarr2_len,tarr2_len;
   register int i;
   int state_var;
   int index;
   int old_triple[ARR_NUM][4];
   int new_triple[ARR_NUM][4];
   int old_len=0;
   int new_len=0;
   int first,second,third;

   zero_neighbours(sarr1,&sarr1_len);
   zero_neighbours(tarr1,&tarr1_len);
   zero_neighbours(sarr2,&sarr2_len);
   zero_neighbours(tarr2,&tarr2_len);

   fill_neighbours(N,s1,t1,1,sarr1,&sarr1_len);
   fill_neighbours(N,t1,s1,0,tarr1,&tarr1_len);
   fill_neighbours(N,s2,t2,1,sarr2,&sarr2_len);
   fill_neighbours(N,t2,s2,0,tarr2,&tarr2_len);

   for (i=0;i<14;i++)
   {
		sub[i]=0;
		add[i]=0;
   }

   for (i=0;i<ARR_NUM;i++)
   {
		old_triple[i][0]=0;old_triple[i][1]=0;old_triple[i][2]=0;old_triple[i][3]=0;
		new_triple[i][0]=0;new_triple[i][1]=0;new_triple[i][2]=0;new_triple[i][3]=0;
   }
   old_len=0;
   new_len=0;

//	fprintf(GNRL_ST.mat_metrop_fp,"starting old set\n");

   /* update sub with the circuits of sarr1,s1,t1 */
   for (i=0;i<sarr1_len;i++)
   {
		refine(old_triple,&old_len,sarr1[i],s1,t1);
   }

   /* update sub with the circuits of s1,t1,tarr1 */
   for (i=0;i<tarr1_len;i++)
   {
		refine(old_triple,&old_len,s1,t1,tarr1[i]);
   }

	/* update sub with the circuits of sarr2,s2,t2 */
   for (i=0;i<sarr2_len;i++)
   {
		refine(old_triple,&old_len,sarr2[i],s2,t2);
   }


   /* update sub with the circuits of s2,t2,tarr2 */
   for (i=0;i<tarr2_len;i++)
   {
		refine(old_triple,&old_len,s2,t2,tarr2[i]);
   }


	/********* New bug fix - when a node connects to both s1 and t2 **********/
   /* update sub with the circuits of sarr1,s1,t2 */
   for (i=0;i<sarr1_len;i++)
   {
		refine(old_triple,&old_len,sarr1[i],s1,t2);
   }

   /* update sub with the circuits of s1,t2,tarr2 */
   for (i=0;i<tarr2_len;i++)
   {
		refine(old_triple,&old_len,s1,t2,tarr2[i]);
   }

	/* update sub with the circuits of sarr2,s2,t1 */
   for (i=0;i<sarr2_len;i++)
   {
		refine(old_triple,&old_len,sarr2[i],s2,t1);
   }


   /* update sub with the circuits of s2,t1,tarr1 */
   for (i=0;i<tarr1_len;i++)
   {
		refine(old_triple,&old_len,s2,t1,tarr1[i]);
   }



   for (i=0;i<old_len;i++)
   {
	   first=old_triple[i][1];
	   second=old_triple[i][2];
	   third=old_triple[i][3];
	   state_var=32*is_edge(N,first,second)+16*is_edge(N,second,third)+8*is_edge(N,third,first)
		+4*is_edge(N,second,first)+2*is_edge(N,third,second)+is_edge(N,first,third);
		index=elm(state_var);
		if ((index>0)&&(index<=13))
			sub[index]=sub[index]+1;
		state_var=32*is_edge2(N,first,second,s1,t1,s2,t2)+16*is_edge2(N,second,third,s1,t1,s2,t2)+8*is_edge2(N,third,first,s1,t1,s2,t2)
		+4*is_edge2(N,second,first,s1,t1,s2,t2)+2*is_edge2(N,third,second,s1,t1,s2,t2)+is_edge2(N,first,third,s1,t1,s2,t2);
		index=elm(state_var);
		if ((index>0)&&(index<=13))
			add[index]=add[index]+1;
   }
}

//   computes the current 13 3-node circuits of s1-t1+neighbours
//   and s2-t2+neighbours. Puts results in sub-a 13*1 vector.
//   computes the future 13 3-node circuits of s1-t2 and neighbours
//   and s2-t1+neighbours. Puts results in add-a 13*1 vector.
//   Takes into account that double edges are switched
void fill_13_dbl(Network *N,int s1,int t1,int s2,int t2,int sub[14],int add[14])

{

   int sarr1[ARR_NUM];
   int tarr1[ARR_NUM];
   int sarr2[ARR_NUM];
   int tarr2[ARR_NUM];
   int sarr1_len,tarr1_len,sarr2_len,tarr2_len;
   register int i;
   int state_var;
   int index;
   int old_triple[ARR_NUM][4];
   int new_triple[ARR_NUM][4];
   int old_len=0;
   int new_len=0;
   int first,second,third;

   zero_neighbours(sarr1,&sarr1_len);
   zero_neighbours(tarr1,&tarr1_len);
   zero_neighbours(sarr2,&sarr2_len);
   zero_neighbours(tarr2,&tarr2_len);

   fill_neighbours(N,s1,t1,1,sarr1,&sarr1_len);
   fill_neighbours(N,t1,s1,0,tarr1,&tarr1_len);
   fill_neighbours(N,s2,t2,1,sarr2,&sarr2_len);
   fill_neighbours(N,t2,s2,0,tarr2,&tarr2_len);

   for (i=0;i<14;i++)
   {
		sub[i]=0;
		add[i]=0;
   }

   for (i=0;i<ARR_NUM;i++)
   {
		old_triple[i][0]=0;old_triple[i][1]=0;old_triple[i][2]=0;old_triple[i][3]=0;
		new_triple[i][0]=0;new_triple[i][1]=0;new_triple[i][2]=0;new_triple[i][3]=0;
   }
   old_len=0;
   new_len=0;

//	fprintf(GNRL_ST.mat_metrop_fp,"starting old set\n");

   /* update sub with the circuits of sarr1,s1,t1 */
   for (i=0;i<sarr1_len;i++)
   {
		refine(old_triple,&old_len,sarr1[i],s1,t1);
   }

   /* update sub with the circuits of s1,t1,tarr1 */
   for (i=0;i<tarr1_len;i++)
   {
		refine(old_triple,&old_len,s1,t1,tarr1[i]);
   }

	/* update sub with the circuits of sarr2,s2,t2 */
   for (i=0;i<sarr2_len;i++)
   {
		refine(old_triple,&old_len,sarr2[i],s2,t2);
   }


   /* update sub with the circuits of s2,t2,tarr2 */
   for (i=0;i<tarr2_len;i++)
   {
		refine(old_triple,&old_len,s2,t2,tarr2[i]);
   }


	/********* New bug fix - when a node connects to both s1 and t2 **********/
   /* update sub with the circuits of sarr1,s1,t2 */
   for (i=0;i<sarr1_len;i++)
   {
		refine(old_triple,&old_len,sarr1[i],s1,t2);
   }

   /* update sub with the circuits of s1,t2,tarr2 */
   for (i=0;i<tarr2_len;i++)
   {
		refine(old_triple,&old_len,s1,t2,tarr2[i]);
   }

	/* update sub with the circuits of sarr2,s2,t1 */
   for (i=0;i<sarr2_len;i++)
   {
		refine(old_triple,&old_len,sarr2[i],s2,t1);
   }


   /* update sub with the circuits of s2,t1,tarr1 */
   for (i=0;i<tarr1_len;i++)
   {
		refine(old_triple,&old_len,s2,t1,tarr1[i]);
   }



   for (i=0;i<old_len;i++)
   {
	   first=old_triple[i][1];
	   second=old_triple[i][2];
	   third=old_triple[i][3];
	   state_var=32*is_edge(N,first,second)+16*is_edge(N,second,third)+8*is_edge(N,third,first)
		+4*is_edge(N,second,first)+2*is_edge(N,third,second)+is_edge(N,first,third);
		index=elm(state_var);
		if ((index>0)&&(index<=13))
			sub[index]=sub[index]+1;
		state_var=32*is_edge2_dbl(N,first,second,s1,t1,s2,t2)+16*is_edge2_dbl(N,second,third,s1,t1,s2,t2)+8*is_edge2(N,third,first,s1,t1,s2,t2)
		+4*is_edge2_dbl(N,second,first,s1,t1,s2,t2)+2*is_edge2_dbl(N,third,second,s1,t1,s2,t2)+is_edge2_dbl(N,first,third,s1,t1,s2,t2);
		index=elm(state_var);
		if ((index>0)&&(index<=13))
			add[index]=add[index]+1;
   }

}


/* checks if there is an edge between node1 and node2 */
int is_edge(Network *N,int node1,int node2)
{
	//register int i;

	return(MatGet(N->mat,node1,node2));

}

/* checks if there is an edge between node1 and node2
   if we are at the cross edges s1->t2 or s2->t1 output 1 */
int is_edge2(Network *N,int node1,int node2,int s1,int t1,int s2,int t2)
{
	//register int i;

	if (((node1==s1)&&(node2==t2))||((node1==s2)&&(node2==t1)))
		return(1);
/*	if (((node2==s1)&&(node1==t2))||((node2==s2)&&(node1==t1)))
		return(1);*/
	if (((node1==s1)&&(node2==t1))||((node1==s2)&&(node2==t2)))
		return(0);
	/*if (((node2==s1)&&(node1==t1))||((node2==s2)&&(node1==t2)))
		return(0);*/
	return(MatGet(N->mat,node1,node2));

}

/* checks if there is an edge between node1 and node2
   if we are at the cross edges s1->t2 or s2->t1 output 1
   considers that a switch means switching putting on
   t2->s1 and t1->s2*/
int is_edge2_dbl(Network *N,int node1,int node2,int s1,int t1,int s2,int t2)
{
	//register int i;

	if (((node1==s1)&&(node2==t2))||((node1==s2)&&(node2==t1)))
		return(1);
	if (((node2==s1)&&(node1==t2))||((node2==s2)&&(node1==t1)))
		return(1);
	if (((node1==s1)&&(node2==t1))||((node1==s2)&&(node2==t2)))
		return(0);
	if (((node2==s1)&&(node1==t1))||((node2==s2)&&(node1==t2)))
		return(0);
	return(MatGet(N->mat,node1,node2));

}

int elm(int state)
{


	switch (state)
		{

			case 10 :
			case 20 :
			case 33 :
					return(1);
			case 11 :
			case 22 :
			case 26 :
			case 37 :
			case 41 :
			case 52 :
					return(2);
			case 27 :
			case 45 :
			case 54 :
					return(3);
			case 12 :
			case 17 :
			case 34 :
					return(4);
			case 3 :
			case 5 :
			case 6 :
			case 24 :
			case 40 :
			case 48 :
					return(5);
			case 13 :
			case 19 :
			case 25 :
			case 38 :
			case 44 :
			case 50 :
					return(6);
			case 7  :
			case 56  :
					return(7);
			case 15 :
			case 23 :
			case 39 :
			case 57 :
			case 58 :
			case 60 :
					return(8);
			case 31 :
			case 47 :
			case 55 :
			case 59 :
			case 61 :
			case 62 :
					return(9);
			case 63 :
					return(10);
			case 14 :
			case 21 :
			case 28 :
			case 35 :
			case 42 :
			case 49 :
					return(11);
			case 29 :
			case 46 :
			case 51 :
					return(12);
			case 30 :
			case 43 :
			case 53 :
					return(13);
			default :
					return(0);

		}
}

/* checks if the triple already occurs in triples. If
   not, adds it. */
void refine(int triples[ARR_NUM][4],int *len,int first,int second,int third)
{
		register int i;
		int t1[4],t2[4],t3[4],t4[4],t5[4],t6[4];
		int check_len;
		int write=1;


		//create all possible 6 permutations of 3 elements

		t1[1]=first; t1[2]=second; t1[3]=third;
		t2[1]=second; t2[2]=third;  t2[3]=first;
		t3[1]=third; t3[2]=first; t3[3]=second;
		t4[1]=first; t4[2]=third; t4[3]=second;
		t5[1]=third; t5[2]=second; t5[3]=first;
		t6[1]=second; t6[2]=first; t6[3]=third;

		check_len=*len;

		if ((first==second)||(second==third)||(first==third))
			return;
		for (i=0;i<check_len;i++)
		{
			if ((equal(t1,triples[i],3))||(equal(t2,triples[i],3))||(equal(t3,triples[i],3))
				||(equal(t4,triples[i],3))||(equal(t5,triples[i],3))||(equal(t6,triples[i],3)))
				//write=0;
				return;

		}

		if (write==1)
		{
			triples[(*len)][1]=first;
			triples[(*len)][2]=second;
			triples[(*len)][3]=third;
			(*len)=(*len)+1;
			//fprintf(GNRL_ST.mat_metrop_fp,"%d %d %d\n",first,second,third);
		}
}


/* returnes 1 if the vectors t1 and t2 are equal between indices 1 and len */
int equal(int *t1,int *t2,int len)
{
	register int i;
	int result=1;

	for (i=1;i<=len;i++)
		if (t1[i]!=t2[i])
			result=0;

	return(result);
}

/* shalev's metropolis additions from main */

/* takes the motif results and puts them in a 13vec */
void met_dump_motifs_res(list64 *res, char *network_name,int vec13[14])
{
	Motif *m;
	list64_item *l_res_id;
	int id,i,index;
	int nadav2elem[239];
	nadav2elem[6]=1;nadav2elem[14]=2;nadav2elem[78]=3;nadav2elem[36]=4;
	nadav2elem[12]=5;nadav2elem[74]=6;nadav2elem[98]=7;nadav2elem[102]=8;nadav2elem[110]=9;
	nadav2elem[238]=10;nadav2elem[38]=11;nadav2elem[108]=12;nadav2elem[46]=13;

	for (i=1;i<=13;i++)
		vec13[i]=0;

	//printf("motifs results for %s :\n", network_name);
	for(l_res_id=list64_get_next(res,NULL); l_res_id!=NULL;
		l_res_id=list64_get_next(res,l_res_id))
		{
			//id=l_res_id->val;
			m=(Motif *)l_res_id->p;
			id=(int)m->id;
			index=nadav2elem[id];
			//vec13[index]=*(int*)l_res_id->p;
			vec13[index]=(int)m->count;
		}
}


void
met_join_subgraphs_res(list64 **res_p, int mtf_sz,int rand_net_indx)
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
	Member *mtf_members_arr, *new_members_arr;
	//double total_prob_counts=0;
	double total_motifs_counts=0;
	list64_init(&final_res);

	//debug
	//for(l_res_id=list_get_next(res,NULL); l_res_id !=NULL; l_res_id=list_get_next(res,l_res_id))
	//	printf("id : %d  occurs : %d \n", l_res_id->val, *(int*)l_res_id->p);


	//pass through res list and for each id look for the isomorphism ids
	while((l_res_id=list64_get_next(res,NULL)) != NULL) {
	//for(l_res_id=list_get_next(res,NULL); l_res_id!=NULL;
	//	l_res_id=list_get_next(res,l_res_id)) {
			iso_id_list=calc_mtf_id_iso((int64)l_res_id->val, mtf_sz);
			//find all other iso id remove them from res list
			//and sum all iso id counts
			//insert new item with iso_id (first id in iso_id_list) to final_res
			l_iso_rep=list64_get_next(iso_id_list,NULL);
			rep_id=(int64)l_iso_rep->val;
			rep_mtf=(Motif*)calloc(1,sizeof(Motif));
			

			rep_mtf->id = rep_id;
			list_init(&rep_mtf->members);
			list_init(&rep_mtf->all_members);

			//pass through iso_id_list
			for(l_iso=list64_get_next(iso_id_list,NULL); l_iso!=NULL;
				l_iso=list64_get_next(iso_id_list,l_iso)) {
					curr_id=(int64)l_iso->val;
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
										new_members_arr=(Member*)calloc(mtf_sz,sizeof(Member));

										//copy it cause later i free all members of l_res motif
										for(i=0;i<mtf_sz;i++){
											new_members_arr[i].node=mtf_members_arr[i].node;
											new_members_arr[i].role=mtf_members_arr[i].role;
										}
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

	//free res list and make res point to final res list
	res_tbl_mem_free_single(res);
	//list64_free_mem(res);
	*res_p=final_res;

}



// compute the real number of motifs
void met_motifs_search_real(Network *N,Res_tbl *met_res_tbl,int vec[14])
{
  // int rc=RC_OK;

	time_measure_start(&GNRL_ST.real_net_time);
	//if(GNRL_ST.run_prob_app==FALSE)
		count_subgraphs(N, 3, &met_res_tbl->real, REAL_NET);
	//else
	//	search_motifs_prob(N, GNRL_ST.mtf_sz, RES_TBL.real, REAL_NET);
	//free_network_mem(N);

	//calc result after isomorphism of ids
	met_join_subgraphs_res(&met_res_tbl->real, 3, 0);
//	fprintf(GNRL_ST.out_fp,"\n   Summary motif results\n");
//	fprintf(GNRL_ST.out_fp,"   =====================\n");

	//Nadav I made the next call comment
	met_dump_motifs_res(met_res_tbl->real,N->name,vec);
	//time_measure_stop(&GNRL_ST.real_net_time);
	//if(GNRL_ST.quiet_mode==FALSE)
	//	dump_time_measure(stdout, "Real network processing runtime was:", &GNRL_ST.real_net_time);

//	return rc;
}



int metrop(double de, double t)
{
	double bolzman,rnd;
	//static long gljdum=1;
	//float frnd;

	bolzman=exp(-de/t);
	rnd=get_rand_double();
	if (de<0.0)
		return(1);
	if ((-(de/t))<-10000.0)
		return(0);
	else
		return(rnd<bolzman);

//	return de < 0.0 || rnd < bolzman;
}

double energy(int target[14],int current[14])
/* computes the energy as sigma(i=1:13){(target[i]-current[i]})/0.5(target[i]+current[i])}*/
{
	register int i;
	double result=0.0;
	//double mone,mechane,sm,df;
	double sm,df;

	for (i=1;i<=13;i++)
	{
		sm=(double)(target[i]+current[i]);
		if (sm!=0)
		{
			df=(double)(abs(target[i]-current[i]));

			result+=df/sm;

			//printf("%lf\t",result);
			//result+=(target[i]-current[i])*(target[i]-current[i]);
		}
	}
	return(result);

}

void update(int target[14],int minus[14],int plus[14],int fwd)
{
	register int i;

	for (i=1;i<=13;i++)
	{
		if (fwd)
			target[i]=target[i]-minus[i]+plus[i];
		else
			target[i]=target[i]+minus[i]-plus[i];

	}
}

void output13(int vec[14],FILE *fp)
{
	register int i;

	//fprintf(fp,"triadic census : ");
	for (i=1;i<=13;i++)
	{
		fprintf(fp,"%d ",vec[i]);
	//	printf("%d ",vec[i]);
	}
	fprintf(fp,"\n");
	//printf("\n");
}


/* perform metropolis iteration exchanging single edges. Output the temperature */
double sin_metrop_switches(Network **RN,int nover,int nlimit,double init_t,int real_vec13[14],int rand_vec13[14])
{

	double t=init_t;
	register int w,k,num_at_t;
	int s1,t1,s2,t2;
	int i,j,l;
	int sub13[14];
	int add13[14];
	double E1,E2,dE;
	int make_change;
	if(GNRL_ST.out_metrop_mat_flag==TRUE)
		printf("\nBeginning single switches\n");

	k=0,w=0,num_at_t=0;
	t=init_t;

	for (i=1;i<14;i++)
	{
		sub13[i]=0;
		add13[i]=0;
	}

	E1=energy(real_vec13,rand_vec13);

	if ((w%1000)==0 && GNRL_ST.out_metrop_mat_flag)
		fprintf(GNRL_ST.mat_metrop_fp,"\nEsin :\n");

	while ((w<nover)&&(E1>GNRL_ST.e_thresh))
	{
		w++;
		num_at_t++;

		if ((k>nlimit)||(num_at_t>nlimit))
		{
			num_at_t=0;
			k=0;
			t *= TFACTR;		// lower temperature
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
			fill_13((*RN),s1,t1,s2,t2,sub13,add13);
			update(rand_vec13,sub13,add13,1);
			E2=energy(real_vec13,rand_vec13);
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
				if(DEBUG_LEVEL>1 && GNRL_ST.out_metrop_mat_flag)
				{
					fprintf(GNRL_ST.mat_metrop_fp,"Changed to  :\n\n");
					for (l=1;l<=(*RN)->e_sin_num;l++)
						fprintf(GNRL_ST.mat_metrop_fp,"%d %d %d\n",(*RN)->e_arr_sin[l].t,(*RN)->e_arr_sin[l].s,1);

					fprintf(GNRL_ST.mat_metrop_fp,"rand_vec : \n");
					output13(rand_vec13,GNRL_ST.mat_metrop_fp);
					fprintf(GNRL_ST.mat_metrop_fp,"- : \n");
					output13(sub13,GNRL_ST.mat_metrop_fp);
					fprintf(GNRL_ST.mat_metrop_fp,"+ : \n");
					output13(add13,GNRL_ST.mat_metrop_fp);
					fprintf(GNRL_ST.mat_metrop_fp,"The network :\n\n");
					for (l=1;l<=(*RN)->e_sin_num;l++)
						fprintf(GNRL_ST.mat_metrop_fp,"%d %d %d\n",(*RN)->e_arr_sin[l].t,(*RN)->e_arr_sin[l].s,1);
					fprintf(GNRL_ST.mat_metrop_fp,"\n\n");
				}
				k++;
			}
			else
			/* change rand_vec13 back */
			{
				update(rand_vec13,sub13,add13,0);

			}
		}
	}
	if(GNRL_ST.out_metrop_mat_flag==TRUE) {
		fprintf(GNRL_ST.mat_metrop_fp,"%lf\n",E1);
		printf("success=%d k=%d E1=%lf T=%lf\n",k,w,E1,t);
	}
	return(t);
}


/* perform metropolis iteration exchanging double edges. Output the temperature */
double dbl_metrop_switches(Network **RN,int nover,int nlimit,double init_t,int real_vec13[14],int rand_vec13[14])
{

	double t=init_t;
	register int w,k,num_at_t;
	int s1,t1,s2,t2;
	int i,j;
	int sub13[14];
	int add13[14];
	double E1,E2,dE;
	int make_change;
	int tries; //number of tries to find proper pair of edges to exchange
	int twin_i=-1,twin_j=-1;
	w=0;
	if ((w%1000)==0 && GNRL_ST.out_metrop_mat_flag)
		fprintf(GNRL_ST.mat_metrop_fp,"\nBeginning double switches\n");
	k=0,w=0,num_at_t=0;
	t=init_t;

	for (i=1;i<14;i++)
	{
		sub13[i]=0;
		add13[i]=0;
	}

	E1=energy(real_vec13,rand_vec13);

	while ((w<nover)&&(E1>GNRL_ST.e_thresh))
	{
		w++;
		num_at_t++;
		if ((k>nlimit)||(num_at_t>nlimit))
		{
			num_at_t=0;
			k=0;
			t *= TFACTR;		// lower temperature
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
			// make fill13_dbl - new function
			fill_13_dbl((*RN),s1,t1,s2,t2,sub13,add13);
			update(rand_vec13,sub13,add13,1);
			E2=energy(real_vec13,rand_vec13);
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
					output13(rand_vec13,GNRL_ST.mat_metrop_fp);
				}
				k++;
			}
			else
			/* change rand_vec13 back */
			{
				update(rand_vec13,sub13,add13,0);
			}
		}
	}
	if(GNRL_ST.out_metrop_mat_flag==TRUE){
		fprintf(GNRL_ST.mat_metrop_fp,"%lf\n",E1);
		printf("success=%d k=%d E1=%lf T=%lf\n",k,w,E1,t);
	}
	return(t);
}


int
gen_rand_network_metrop(Network **RN_p,int real_vec13[14])
{
	int rc=RC_OK;
	int sin_num_of_switchs,dbl_num_of_switchs,num_at_t;
	//int s2,t1,t2;
	int i,j,l;
	int k;  // number of succesful switches
	int w;  // total number of switches
	//int tries; //number of tries to find proper pair of edges to exchange
	//int twin_i,twin_j;
	Network *RN;
	int rand_network_checked=FALSE;
	Res_tbl met_res_tbl;
	//int real_vec13[14];
	int rand_vec13[14];
	//int sub13[14];
	//int add13[14];
	//double E1,E2,dE;
	double E1,E2;
	int	snover,dnover; // maximum # of graph changes at any temperature
	//int nlimit,snlimit,dnlimit; // maximum # of succesful graph changes at any temperature
	int snlimit,dnlimit; // maximum # of succesful graph changes at any temperature
	//int make_change;
	//int numsucc; // how many successful switches already made
	double t;  // the temperature
	int success;
	int dispair=0;
	double dummy;
	int switches_range;
	double *switch_ratio_arr;
	switch_ratio_arr=(double*)calloc(GNRL_ST.rnd_net_num+1,sizeof(double));

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

		gen_rand_network_switches_method_conserve_double_edges(&RN, &dummy);



	/**************************************************************/
	/**************************************************************/
	/**************************************************************/
	/* Part II - metropolis algorithm. Perform switches with probability*/
	/**************************************************************/
	/**************************************************************/
	/**************************************************************/

	/**************************************************************/
		/* metropolis parameters depend on num_of _switches*/


	/* otuput real network to _METROP file */



	init_res_tbl(&met_res_tbl);
	list64_init(&met_res_tbl.real);
	met_motifs_search_real(RN,&met_res_tbl,rand_vec13);
	list64_free_mem(met_res_tbl.real);

	E1=energy(real_vec13,rand_vec13);

	w=0;
	if(GNRL_ST.out_metrop_mat_flag) {
		fprintf(GNRL_ST.mat_metrop_fp,"\nreal network\n");
		fprintf(GNRL_ST.mat_metrop_fp,"E=%lf T=%lf\n",E1,t);
		fprintf(GNRL_ST.mat_metrop_fp,"Real triadic census :\n ");
		output13(real_vec13,GNRL_ST.mat_metrop_fp);
		fprintf(GNRL_ST.mat_metrop_fp,"random network\n");
		fprintf(GNRL_ST.mat_metrop_fp,"%lf,%lf \n",E1,t);
		fprintf(GNRL_ST.mat_metrop_fp,"Initial random triadic census :\n ");
		output13(rand_vec13,GNRL_ST.mat_metrop_fp);
	}

/*	nover=sin_num_of_switchs*GNRL_ST.iteration_factor;
	nlimit=nover/10;
	k=0;
	num_at_t=0;
	sin_metrop_switches(&RN,nover,nlimit,GNRL_ST.t_init,real_vec13,rand_vec13);
*/
	if (dbl_num_of_switchs==0)
	{
		snover=sin_num_of_switchs*(int)GNRL_ST.iteration_factor;
		snlimit=snover/10;
		t=sin_metrop_switches(&RN,snover,snlimit,GNRL_ST.t_init,real_vec13,rand_vec13);
	}
	else
	{
		dnover=dbl_num_of_switchs*(int)GNRL_ST.iteration_factor/10;
		dnlimit=dnover/10;
		if (dnover<10)
		{

			dnover=dbl_num_of_switchs;
			dnlimit=dnover;
		}
		t=dbl_metrop_switches(&RN,dnover,dnlimit,GNRL_ST.t_init,real_vec13,rand_vec13);

		snover=sin_num_of_switchs*(int)GNRL_ST.iteration_factor/10;
		snlimit=snover/10;
		k=0;
		num_at_t=0;
		t=sin_metrop_switches(&RN,snover,snlimit,t,real_vec13,rand_vec13);

		if(GNRL_ST.out_metrop_mat_flag) {
			fprintf(GNRL_ST.mat_metrop_fp,"End E\n");
			printf("End E\n");
		}
		for (i=1;i<10;i++)
		{
			E2=energy(real_vec13,rand_vec13);
			if (E2>0)
			{
				t=dbl_metrop_switches(&RN,dnover,dnlimit,t,real_vec13,rand_vec13);
				t=sin_metrop_switches(&RN,snover,snlimit,t,real_vec13,rand_vec13);
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
		output13(real_vec13,GNRL_ST.mat_metrop_fp);
		fprintf(GNRL_ST.mat_metrop_fp,"Final random triadic census :\n ");
		output13(rand_vec13,GNRL_ST.mat_metrop_fp);
	}

	*RN_p=RN;
	return RC_OK;
}
