/************************************************************************
*
*  File name: motif_ids.c
*
*  Description: Motif ids interface and functions
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "common.h"


/*************************** Global variables ****************************/
static Matrix *row_p_m, *col_p_m, *id_m, *res_m, *tmp_m;


/******************************* Externs *********************************/

extern int DEBUG_LEVEL;
extern Gnrl_st GNRL_ST;
extern int p2[2][2];
extern int p3[6][3];
extern int p4[24][4];
extern int p5[120][5];
extern int p6[720][6];
extern int p7[5040][7];
extern int p8[40320][8];

/******************************* Functions *******************************/
int *
get_perm(int size,int index){
	switch (size){
	case 2:
		return p2[index];
		break;
	case 3:
		return p3[index];
		break;
	case 4:
		return p4[index];
		break;
	case 5:
		return p5[index];
		break;
	case 6:
		return p6[index];
		break;
	case 7:
		return p7[index];
		break;
	case 8:
		return p8[index];
		break;
	default:
		return NULL;
	}
}

int
gen_subset_matrix(Network *N, Network **SN, list *vrtx_set, int mtf_sz)
{
	int i,j,rc=RC_OK;
	list_item *row,*col;
	(*SN)=(Network*)calloc(1,sizeof(Network));
	(*SN)->vertices_num=mtf_sz;
	rc |= MatInit(&((*SN)->mat),(*SN)->vertices_num,FULL);
	if (rc == RC_ERR) {
		printf("Malloc failed!!!\n");
		return RC_ERR;
	}

	(*SN)->edges_num = 0;
	(*SN)->e_arr=NULL;
	(*SN)->name=NULL;

	row=vrtx_set->l;
	col=vrtx_set->l;
	//make sure all matrix is init to zeros
	for (i=1; i<=(*SN)->vertices_num; i++) {
		for (j=1; j<=(*SN)->vertices_num; j++) {
			MatAsgn((*SN)->mat,i,j,MatGet(N->mat,row->val,col->val));
			col=col->next;
		}
		row=row->next;
		col=vrtx_set->l;
	}
	return RC_OK;
}


void
fill_mat_id(Matrix *M, int64 id)
{
	int i;
	int mat_cells_num=M->size*M->size;
	for(i=0;i<mat_cells_num;i++,id/=2)
		M->m[i]=((char)(id & 1));
}

void
static fill_row_perm_mat(Matrix *M, int *p_arr)
{
	int i;
	for (i=0;i<M->size;i++)
		MTRX(M,i+1,p_arr[i]+1)=1;
}
void
static fill_col_perm_mat(Matrix *M, int *p_arr)
{
	int i;
	for (i=0;i<M->size;i++)
		MTRX(M,p_arr[i]+1,i+1)=1;
}


//multiply square matrices
void
static multiply_mat(Matrix *L_M, Matrix *R_M, Matrix *RES_M)
{
	int i,j,k;

	for (i=1;i<=L_M->size;i++)
		for(j=1;j<=L_M->size;j++) {
			for(k=1;k<=L_M->size;k++)
				MTRX(RES_M,i,j)+=(char)MTRX(L_M,i,k)*(char)MTRX(R_M,k,j);
		}
}





int64
static get_motif_id(Matrix *M)
{
	int i;
	int mat_cells_num=M->size*M->size;
	int64 id=0;

	//id is the binary of the sub-matrix seen as one long array
	for(i=0;i<mat_cells_num;i++)
		id+=(int64)pow(2,i)*(int64)(M->m[i]);
	return id;
}


list64*
calc_id_iso_2(int64 id)
{
	int i;
	list64 *iso_list;
	int64 iso_id;
	int perm_num=2;
	int *perm_idx;

	list64_init(&iso_list);
	for(i=0; i<perm_num; i++) {
		fill_row_perm_mat(row_p_m, (int*)p2[i]);
		//dump_matrix(row_p_m, "row_perm");
		fill_col_perm_mat(col_p_m, (int*)p2[i]);
		//dump_matrix(col_p_m, "col_perm");
		fill_mat_id(id_m,id);
		//dump_matrix(id_m, "id_m");
		multiply_mat(row_p_m, id_m, tmp_m);
		//dump_matrix(tmp_m, "tmp_m");
		multiply_mat(tmp_m, col_p_m, res_m);
		//dump_matrix(res_m, "res_m");
		//calc id from matrix

		iso_id=get_motif_id(res_m);
		//insert new id to iso list (if not already there)
		if(list64_get(iso_list,iso_id)==NULL){
			perm_idx=(int*)calloc(1,sizeof(int));
			*perm_idx=i;
			list64_insert(iso_list,iso_id,(void*)perm_idx);
		}
		clear_matrix(row_p_m);
		clear_matrix(col_p_m);
		clear_matrix(tmp_m);
		clear_matrix(res_m);
	}
	return iso_list;
}

list64*
calc_id_iso_3(int64 id)
{
	int i;
	list64 *iso_list;
	int64 iso_id;
	int perm_num=6;
	int *perm_idx;

	list64_init(&iso_list);
	for(i=0; i<perm_num; i++) {
		fill_row_perm_mat(row_p_m, (int*)p3[i]);
		//dump_matrix(row_p_m, "row_perm");
		fill_col_perm_mat(col_p_m, (int*)p3[i]);
		//dump_matrix(col_p_m, "col_perm");
		fill_mat_id(id_m,id);
		//dump_matrix(id_m, "id_m");
		multiply_mat(row_p_m, id_m, tmp_m);
		//dump_matrix(tmp_m, "tmp_m");
		multiply_mat(tmp_m, col_p_m, res_m);
		//dump_matrix(res_m, "res_m");

		//calc id from matrix
		iso_id=get_motif_id(res_m);

		//insert new id to iso list (if not already there)
		if(list64_get(iso_list,iso_id)==NULL){
			perm_idx=(int*)calloc(1,sizeof(int));
			*perm_idx=i;
			list64_insert(iso_list,iso_id,(void*)perm_idx);
		}
		clear_matrix(row_p_m);
		clear_matrix(col_p_m);
		clear_matrix(tmp_m);
		clear_matrix(res_m);
	}
	return iso_list;
}

list64*
calc_id_iso_4(int64 id)
{
	int i;
	list64 *iso_list;
	int64 iso_id;
	int perm_num=24;
	int *perm_idx;

	list64_init(&iso_list);
	for(i=0; i<perm_num; i++) {
		fill_row_perm_mat(row_p_m, (int*)p4[i]);
		fill_col_perm_mat(col_p_m, (int*)p4[i]);
		fill_mat_id(id_m,id);
		multiply_mat(row_p_m, id_m, tmp_m);
		multiply_mat(tmp_m, col_p_m, res_m);
		//calc id from matrix
		iso_id=get_motif_id(res_m);
		//insert new id to iso list (if not already there)
		if(list64_get(iso_list,iso_id)==NULL){
			perm_idx=(int*)calloc(1,sizeof(int));
			*perm_idx=i;
			list64_insert(iso_list,iso_id,(void*)perm_idx);
		}
		clear_matrix(row_p_m);
		clear_matrix(col_p_m);
		clear_matrix(tmp_m);
		clear_matrix(res_m);
	}
	return iso_list;
}

list64*
calc_id_iso_5(int64 id)
{
	int i;
	list64 *iso_list;
	int64 iso_id;
	int perm_num=120;
	int *perm_idx;

	list64_init(&iso_list);
	for(i=0; i<perm_num; i++) {
		fill_row_perm_mat(row_p_m, (int*)p5[i]);
		fill_col_perm_mat(col_p_m, (int*)p5[i]);
		fill_mat_id(id_m,id);
		multiply_mat(row_p_m, id_m, tmp_m);
		multiply_mat(tmp_m, col_p_m, res_m);
		//calc id from matrix
		iso_id=get_motif_id(res_m);
		//insert new id to iso list (if not already there)
		if(list64_get(iso_list,iso_id)==NULL){
			perm_idx=(int*)calloc(1,sizeof(int));
			*perm_idx=i;
			list64_insert(iso_list,iso_id,(void*)perm_idx);
		}
		clear_matrix(row_p_m);
		clear_matrix(col_p_m);
		clear_matrix(tmp_m);
		clear_matrix(res_m);
	}
	return iso_list;
}

list64*
calc_id_iso_6(int64 id)
{
	int i;
	list64 *iso_list;
	int64 iso_id;
	int perm_num=720;
	int *perm_idx;

	list64_init(&iso_list);
	for(i=0; i<perm_num; i++) {
		fill_row_perm_mat(row_p_m, (int*)p6[i]);
		fill_col_perm_mat(col_p_m, (int*)p6[i]);
		fill_mat_id(id_m,id);
		multiply_mat(row_p_m, id_m, tmp_m);
		multiply_mat(tmp_m, col_p_m, res_m);
		//calc id from matrix
		iso_id=get_motif_id(res_m);
		//insert new id to iso list (if not already there)
		if(list64_get(iso_list,iso_id)==NULL){
			perm_idx=(int*)calloc(1,sizeof(int));
			*perm_idx=i;
			list64_insert(iso_list,iso_id,(void*)perm_idx);
		}
		clear_matrix(row_p_m);
		clear_matrix(col_p_m);
		clear_matrix(tmp_m);
		clear_matrix(res_m);
	}
	return iso_list;
}

list64*
calc_id_iso_7(int64 id)
{
	int i;
	list64 *iso_list;
	int64 iso_id;
	int perm_num=5040;
	int *perm_idx;

	list64_init(&iso_list);
	for(i=0; i<perm_num; i++) {
		fill_row_perm_mat(row_p_m, (int*)p7[i]);
		fill_col_perm_mat(col_p_m, (int*)p7[i]);
		fill_mat_id(id_m,id);
		multiply_mat(row_p_m, id_m, tmp_m);
		multiply_mat(tmp_m, col_p_m, res_m);
		//calc id from matrix
		iso_id=get_motif_id(res_m);
		//insert new id to iso list (if not already there)
		if(list64_get(iso_list,iso_id)==NULL){
			perm_idx=(int*)calloc(1,sizeof(int));
			*perm_idx=i;
			list64_insert(iso_list,iso_id,(void*)perm_idx);
		}
		clear_matrix(row_p_m);
		clear_matrix(col_p_m);
		clear_matrix(tmp_m);
		clear_matrix(res_m);
	}
	return iso_list;
}

list64*
calc_id_iso_8(int64 id)
{
	int i;
	list64 *iso_list;
	int64 iso_id;
	int perm_num=40320;
	int *perm_idx;

	list64_init(&iso_list);
	for(i=0; i<perm_num; i++) {
		fill_row_perm_mat(row_p_m, (int*)p8[i]);
		fill_col_perm_mat(col_p_m, (int*)p8[i]);
		fill_mat_id(id_m,id);
		multiply_mat(row_p_m, id_m, tmp_m);
		multiply_mat(tmp_m, col_p_m, res_m);
		//calc id from matrix
		iso_id=get_motif_id(res_m);
		//insert new id to iso list (if not already there)
		if(list64_get(iso_list,iso_id)==NULL){
			perm_idx=(int*)calloc(1,sizeof(int));
			*perm_idx=i;
			list64_insert(iso_list,iso_id,(void*)perm_idx);
		}
		clear_matrix(row_p_m);
		clear_matrix(col_p_m);
		clear_matrix(tmp_m);
		clear_matrix(res_m);
	}
	return iso_list;
}


//get an id as an input and output all identical
//motif ids
//return list of all id isomorfisms
list64*
calc_mtf_id_iso(int64 id, int mtf_sz)
{
	list64 *iso_list=NULL;
	
	row_p_m=init_matrix(mtf_sz);
	col_p_m=init_matrix(mtf_sz);
	id_m=init_matrix(mtf_sz);
	res_m=init_matrix(mtf_sz);
	tmp_m=init_matrix(mtf_sz);
	
	switch (mtf_sz) {
	case 2:
		iso_list=calc_id_iso_2(id);
		break;
	case 3:
		iso_list=calc_id_iso_3(id);
		break;
	case 4:
		iso_list=calc_id_iso_4(id);
		break;
	case 5:
		iso_list=calc_id_iso_5(id);
		break;
	case 6:
		iso_list=calc_id_iso_6(id);
		break;
	case 7:
		iso_list=calc_id_iso_7(id);
		break;
	case 8:
		iso_list=calc_id_iso_8(id);
		break;
	default:
		break;
	}
	free_matrix(row_p_m);
	free(row_p_m);
	free_matrix(col_p_m);
	free(col_p_m);
	free_matrix(id_m);
	free(id_m);
	free_matrix(res_m);
	free(res_m);
	free_matrix(tmp_m);
	free(tmp_m);

	return iso_list;
}


int64
get_rep_mtf_id_3(int64 mtf_id)
{
	switch ((int)mtf_id) {
	case 6:
	case 40:
	case 192:
		return 6;
	case 12:
	case 34:
	case 66:
	case 96:
	case 132:
	case 136:
		return 12;
	case 14:
	case 42:
	case 70:
	case 168:
	case 196:
	case 224:
		return 14;
	case 36:
	case 72:
	case 130:
		return 36;
	case 38:
	case 44:
	case 104:
	case 134:
	case 194:
	case 200:
		return 38;
	case 46:
	case 198:
	case 232:
		return 46;
	case 74:
	case 76:
	case 100:
	case 138:
	case 162:
	case 164:
		return 74;
	case 78:
	case 170:
	case 228:
		return 78;
	case 98:
	case 140:
		return 98;
	case 102:
	case 106:
	case 142:
	case 172:
	case 204:
	case 226:
		return 102;
	case 108:
	case 166:
	case 202:
		return 108;
	case 110:
	case 174:
	case 206:
	case 230:
	case 234:
	case 236:
		return 110;
	case 238:
		return 238;
	default:
		return 0;
	}
}



