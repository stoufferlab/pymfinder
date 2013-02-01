/************************************************************************
*
*  File name: mat.h
*
*  Description: Header file
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#ifndef __MAT_H
#define __MAT_H


#include <stdio.h>
#include "list.h"


/******************* Definitions ************************/


#define FULL          0
#define SPARSE        1

//maximum size of full matrix representation
#define MAT_MAX_SIZE  9


/******************* Structures *************************/


typedef struct {
	int size; //size of matrix (num of rows/collumns)
	char *m;  //full matrix - 2 dim array of chars
}Matrix;

typedef struct{
	list *to;
	list *from;
	int  self_edge; /* 1 if there is a self edge in the edge_list, 0 - otherwise */
}Edges_lists;

typedef struct {
	int size; //size of matrix (num of rows/collumns)
	Edges_lists *m; //sparse matrix - array of lists of edges
}SparseMatrix;

typedef struct {
	int type;  //FULL or SPARSE
	union {
	Matrix *full; //full matrix (values are chars only)
	SparseMatrix *spr; //sparse matrix (values are integers)
	};
}Mat;



/******************* Macros *****************************/


//M should be of type Matrix
#define MTRX(M,i,j)      ( *(char*)(M->m+((i-1)*(M->size))+(j-1)) )




/******************* Prototypes *************************/

//full matrix function
Matrix*
init_matrix(int size);
void
clear_matrix(Matrix *M);
void
free_matrix(Matrix *M);
void
dump_matrix(FILE *FP, Matrix *M, char *name);


//general matrix functions
int
MatGet(Mat *mat,int i, int j);

//asign valuse val to mat[i,j]
void
MatAsgn(Mat *mat, int i, int j, int val);

int
MatInit(Mat **mat_p,int size,int type);

void
MatFree(Mat *mat);






//MACROS

//#define MAT(N,i,j)      ( *(char*)(N->m+((i-1)*(N->vertices_num))+(j-1)) )


#endif
