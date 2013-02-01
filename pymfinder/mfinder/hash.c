/************************************************************************
*
*  File name: hash.c
*
*  Description: Hash functions
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "list.h"
#include "hash.h"


/*************************** Global variables ****************************/


/******************************* Externs *********************************/
extern int DEBUG_LEVEL;
extern Gnrl_st GNRL_ST;



//qsort compare
int compare( const void *arg1, const void *arg2 )
{
   //compare ints
   return ( *(int*)arg1 > *(int*)arg2 );
}



void
hash_init(Hash*hash,int size,int bckt_sz)
{
	hash->size=size;
	hash->bckt_sz=bckt_sz;
	hash->h=(list**)calloc((size),sizeof(list*));
}


Hash*
multi_hash_init(int hash_tbl_num)
{
	int i;
	Hash* multi_hash;

	multi_hash = (Hash*)calloc(hash_tbl_num+1,sizeof(Hash));
	for(i=1;i<=hash_tbl_num;i++)
		hash_init(&multi_hash[i],HASH_SIZE,i);
	return multi_hash;
}

//hush function is a polynom value mod 'hash size'
//the polynom is of a base 'HASH_KEY_BASE'
//when each item in the bucket is the coefficient of the right power of the base.
//
int
hash_get_key(Hash *hash,int *val)
{
	double poly=0;
	int i;
	int key;

	if(hash->bckt_sz==1)
		poly=(double)*val;
	else
		for(i=0; i<hash->bckt_sz; i++)
			poly = poly + (double)(val[i]*pow((double)HASH_KEY_BASE,i));
	key=(int)fmod(poly,(double)hash->size);
	return key;
}



void*
hash_get(Hash *hash, int val)
{
	int key;
	list_item *tmp;

	key=hash_get_key(hash,&val);
	if( (key<0) || (key>=hash->size) ) {
		printf("\nError: Hash function Error\n");
		at_exit(-1);
	}
	if (hash->h[key] == NULL)
		return NULL;

	tmp = list_get(hash->h[key],val);

	if(tmp!=NULL)
		return tmp->p;
	else
		return NULL;
}

void
hash_insert(Hash *hash, void*item, int*val)
{
	int key;

	if(val==NULL)
		return;
	if(hash->bckt_sz>1)
		//sort val array
		qsort( (void *)val, (size_t)hash->bckt_sz, sizeof(int), compare );

	key=hash_get_key(hash,val);
	if(hash->h[key]==NULL)
		list_init(&hash->h[key]);
	list_insert(hash->h[key],(int)*val,item);
}

void
hash_free_mem(Hash *hash)
{
	int i;
	for(i=0;i<hash->size;i++){
		if(hash->h[i]!=NULL) {
			list_free_mem(hash->h[i]);
			//free((void*)hash->h[i]);
		}
	}
	free((void*)hash->h);
}

void
multi_hash_free_mem(Hash *multi_hash, int hash_tbl_num)
{
	int i;

	for(i=1; i<hash_tbl_num+1; i++)
		hash_free_mem(&multi_hash[i]);
	free((void*)multi_hash);
}



