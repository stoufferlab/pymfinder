/************************************************************************
*
*  File name: mat.c
*
*  Description: matrix operations functions
*  
*  Copyright © 2002-2004 Weizmann Institute of Science, 
*			   76100 Rehovot Israel, All rights reserved 
*
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mat.h"


/*************************** Global variables ****************************/


/******************************* Externs *********************************/


extern int DEBUG_LEVEL;



/******************************* Functions *******************************/




Matrix*
init_matrix(int size)
{
	Matrix *mat;
	mat=(Matrix*)calloc(1,sizeof(Matrix));
	//if (DEBUG_LEVEL >=9)
				//fprintf(GNRL_ST.log_fp,"function:init_matrix : mat allocated at %x\n",mat);

	mat->size=size;
	mat->m=(char*)calloc((size+1)*(size+1), sizeof(char));
	//if (DEBUG_LEVEL >=9)
				//fprintf(GNRL_ST.log_fp,"function:init_matrix : mat->m allocated at %x\n",mat->m);
	return mat;
}

SparseMatrix*
init_spr_matrix(int size)
{
	SparseMatrix *mat;
	mat=(SparseMatrix*)calloc(1,sizeof(SparseMatrix));
	mat->size=size;
	mat->m=(Edges_lists*)calloc((unsigned int)size+1,sizeof(Edges_lists));
	return mat;
}

void
clear_matrix(Matrix *M)
{
	int i;
	int mat_cells_num=M->size*M->size;

	for(i=0;i<mat_cells_num;i++)
		M->m[i]=0;
}

void
free_matrix(Matrix *M)
{
	if (M != NULL) {
		free(M->m); 
		//free((void*)M);
	}
}

void
free_spr_matrix(SparseMatrix *M)
{
	int i;

	if(M != NULL) {
		for(i=0;i<=M->size;i++) {
			if(M->m[i].to != NULL)
				list_free_mem(M->m[i].to);
			if(M->m[i].from != NULL)
				list_free_mem(M->m[i].from);
		}
		free(M->m);
	}
}
	


void
dump_matrix(FILE *fp,Matrix *M, char *name)
{
	int i,j;
	if(strcmp(name,""))
		fprintf(fp,"name :%s\n",name);

	for (i=1; i<=M->size; i++) {
		fprintf(fp,"\n");
		for (j=1; j<=M->size; j++)
			fprintf (fp,"%d ",MTRX(M,i,j));
	}
	fprintf(fp,"\n");
}

void
dump_spr_matrix(FILE *fp,Mat *M, char *name)
{
	int i,j;
	if(strcmp(name,""))
		fprintf(fp,"name :%s\n",name);

	for (i=1; i<=M->spr->size; i++) {
		fprintf(fp,"\n");
		for (j=1; j<=M->spr->size; j++)
			fprintf (fp,"%d ",MatGet(M,i,j));
	}
	fprintf(fp,"\n");
}



//General matrices functions - full or sparse

//returns value of mat[i,j]
int
MatGet(Mat *mat,int i, int j)
{
	list_item *e_l;

	switch (mat->type) {
	case FULL:
		return MTRX(mat->full,i,j);
		break;
	case SPARSE:
		if(i<0 || j<0)
			printf("Error in matget\n");
		if(i == j)
		{
			return (mat->spr->m[i].self_edge);
		}
		if( (e_l=list_get(mat->spr->m[i].to,j)) != NULL)
			//edge exist
			return *(int*)e_l->p;
		else
			//edge doest not exist return zero
			return 0;
		break;
	default:
		return -1;
		break;
	}
}

//asign values val to mat[i,j]
void
MatAsgn(Mat *mat, int i, int j, int val)
{
	list_item *e_l;
	int *val_p;

	switch (mat->type) {
	case FULL:
		MTRX(mat->full,i,j)=(char)val;
		break;
	case SPARSE:
		if(i == j)
		{
			mat->spr->m[i].self_edge = 1; //self edge
			break;
		}
		//update the "to" lists
		if( mat->spr->m[i].to==NULL){
			//no edges from this vertex
			list_init(&mat->spr->m[i].to);
			val_p=(int*)calloc(1,sizeof(int));
			*val_p=val;
			list_insert(mat->spr->m[i].to,j,val_p);
		}else if( (e_l=list_get(mat->spr->m[i].to,j)) == NULL) {
			//edges from this vertex exist but not this one
			val_p=(int*)calloc(1,sizeof(int));
			*val_p=val;
			list_insert(mat->spr->m[i].to,j,val_p);
		}else{
			//edge exist but need to assign val
			if (val==0)
				//need to delete this edge from list
				list_delete(mat->spr->m[i].to,j);
			else
				//different value to be assigned
				*(int*)e_l->p=val;
		}
		//update the "from" lists
		if( mat->spr->m[j].from==NULL){
			//no edges from this vertex
			list_init(&mat->spr->m[j].from);
			val_p=(int*)calloc(1,sizeof(int));
			*val_p=val;
			list_insert(mat->spr->m[j].from,i,val_p);
		}else if( (e_l=list_get(mat->spr->m[j].from,i)) == NULL) {
			//edges from this vertex exist but not this one
			val_p=(int*)calloc(1,sizeof(int));
			*val_p=val;
			list_insert(mat->spr->m[j].from,i,val_p);
		}else{
			//edge exist but need to assign val
			if (val==0)
				//need to delete this edge from list
				list_delete(mat->spr->m[j].from,i);
			else
				//different value to be assigned
				*(int*)e_l->p=val;
		}
		break;
	default:
		break;
	}
}

//arguments:
//type = FULL/SPARSE
int
MatInit(Mat **mat_p, int size, int type)
{
	Mat *mat;
	mat=(Mat*)calloc(1,sizeof(Mat));
	//if (DEBUG_LEVEL >=9)
				//fprintf(GNRL_ST.log_fp,"function MatInit: mat allocated at %x\n",mat);
	
	//if number of vertices not bigger then MAT_MAX_SIZE
	//then use full matrix data structure
	//else use sparse matrix data structure
	//if(size <= MAT_MAX_SIZE)
		//mat->type=FULL;
	//else
		//mat->type=SPARSE;
	mat->type=type;

	switch(mat->type){
	case FULL:
		mat->full=init_matrix(size);
		break;
	case SPARSE:
		mat->spr=init_spr_matrix(size);
		break;
	default:
		return -1;
		break;
	}
	*mat_p=mat;
	return 0;
}

void
MatFree(Mat *mat)
{
	switch(mat->type){
	case FULL:
		free_matrix(mat->full);
		free((void*)mat->full);
		break;
	case SPARSE:
		free_spr_matrix(mat->spr);
		free((void*)mat->spr);
		break;
	default:
		break;
	}
}




