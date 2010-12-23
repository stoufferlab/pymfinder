/************************************************************************
*
*  File name: list.c
*
*  Description: list operations functions
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "list.h"

/*************************** Global variables ****************************/


/******************************* Externs *********************************/


/******************************* Functions *******************************/
extern int DEBUG_LEVEL;

void
list_init(list **L)
{
	list *l_new;

	l_new = (list*)calloc(1,sizeof(list));
	l_new->size=0;
	l_new->l=NULL;

	(*L)=l_new;
}


void
list64_init(list64 **L)
{
	list64 *l_new;

	l_new = (list64*)calloc(1,sizeof(list64));
	l_new->size=0;
	l_new->l=NULL;

	(*L)=l_new;
}

void
list_insert(list *L, int val, void *p)
{
	list_item *item,*tmp;
	int pos_found=0;
	item = (list_item*)calloc(1,sizeof(list_item));
	item->val = val;
	item->p=p;
	item->next=NULL;

	//if empty list
	if (L->l == NULL) {
		L->l = item;
	} else {
		tmp = L->l;
		//if smaller then first item
		if (item->val <= L->l->val) {
			item->next=L->l;
			L->l=item;
		} else {
			while( tmp->next!=NULL && !pos_found){
			if( tmp->next->val <= item->val )
				tmp=tmp->next;
			else
				pos_found=1;
			}
			item->next=tmp->next;
			tmp->next=item;
		}
	}
	L->size++;
}

void
list64_insert(list64 *L, int64 val, void *p)
{
	list64_item *item,*tmp;
	int pos_found=0;
	item = (list64_item*)calloc(1,sizeof(list64_item));
	item->val = val;
	item->p=p;
	item->next=NULL;

	//if empty list
	if (L->l == NULL) {
		L->l = item;
	} else {
		tmp = L->l;
		//if smaller then first item
		if (item->val <= L->l->val) {
			item->next=L->l;
			L->l=item;
		} else {
			while( tmp->next!=NULL && !pos_found){
			if( tmp->next->val <= item->val )
				tmp=tmp->next;
			else
				pos_found=1;
			}
			item->next=tmp->next;
			tmp->next=item;
		}
	}
	L->size++;
}

int
list_delete(list *L, int val)
{
	list_item *t1,*t2;

	//if empty list
	if (L->l == NULL) {
		return -1;
	} else {
		t1 = L->l;
		t2 = t1->next;
		//if first item
		if (t1->val == val) {
			L->l = t1->next;
			L->size--;
			if (t1->p != NULL)
				free(t1->p);
			free(t1);
			return 0;
		}
		if (t2 != NULL) {
			while(t2->next != NULL && t2->val!=val) {
				t1=t1->next;
				t2=t2->next;
			}
			if (t2->val == val) {
				t1->next=t2->next;
				L->size--;
				if (t2->p != NULL)
					free(t2->p);
				free(t2);
			} else
				//not in the list
				return -1;
		}
		else return -1;
	}
	return 0;
}

int
list64_delete(list64 *L, int64 val)
{
	list64_item *t1,*t2;

	//if empty list
	if (L->l == NULL) {
		return -1;
	} else {
		t1 = L->l;
		t2 = t1->next;
		//if first item
		if (t1->val == val) {
			L->l = t1->next;
			L->size--;
			if (t1->p != NULL)
				free(t1->p);
			free(t1);
			return 0;
		}
		if (t2 != NULL) {
			while(t2->next != NULL && t2->val!=val) {
				t1=t1->next;
				t2=t2->next;
			}
			if (t2->val == val) {
				t1->next=t2->next;
				L->size--;
				if (t2->p != NULL)
					free(t2->p);
				free(t2);
			} else
				//not in the list
				return -1;
		}
		else return -1;
	}
	return 0;
}

list_item*
list_get(list *L, int val)
{
	list_item *tmp;
	//if list is NULL
	if (L == NULL)
		return NULL;
	//if empty list
	if (L->l == NULL) {
		return NULL;
	} else
		tmp = L->l;

	while(tmp->next != NULL && tmp->val!=val){
		tmp=tmp->next;
	}
	if (tmp->val == val)
		return tmp;
	else
	//not in the list
		return NULL;
}

list64_item*
list64_get(list64 *L, int64 val)
{
	list64_item *tmp;

	//if list is NULL
	if (L == NULL)
		return NULL;
	//if empty list
	if (L->l == NULL) {
		return NULL;
	} else
		tmp = L->l;

	while(tmp->next != NULL && tmp->val!=val){
		tmp=tmp->next;
	}
	if (tmp->val == val)
		return tmp;
	else
	//not in the list
		return NULL;
}


list_item*
list_get_by_indx(list *L, int indx)
{
	list_item *tmp;
	int i;

	//if list is NULL
	if (L == NULL)
		return NULL;
	//if empty list
	if (L->l == NULL) {
		return NULL;
	} else
		tmp = L->l;

	if(indx==0)
		return NULL;

	for(tmp=L->l,i=1;(tmp!=NULL) && (i<indx);tmp=tmp->next,i++);

	if (i==indx)
		return tmp;
	else
		return NULL;
}


list64_item*
list64_get_by_indx(list64 *L, int indx)
{
	list64_item *tmp;
	int i;

	//if list is NULL
	if (L == NULL)
		return NULL;
	//if empty list
	if (L->l == NULL) {
		return NULL;
	} else
		tmp = L->l;

	if(indx==0)
		return NULL;

	for(tmp=L->l,i=1;(tmp!=NULL) && (i<indx);tmp=tmp->next,i++);

	if (i==indx)
		return tmp;
	else
		return NULL;
}



//return next element in list.
// to get first item second parameter should be NULL
//if end of list returns NULL
list_item *
list_get_next(list *L, list_item *item)
{
	if(L==NULL)
		return NULL;
	if(item==NULL)
		return L->l;
	else
		return item->next;
}

//return next element is list.
// to get first item second parameter should be NULL
//if end of list returns NULL
list64_item *
list64_get_next(list64 *L, list64_item *item)
{
	if(L==NULL)
		return NULL;
	if(item==NULL)
		return L->l;
	else
		return item->next;
}



void
list_free_mem(list *L)
{
	list_item *item,*tmp;
	if(L==NULL)
		return;
	for(item=list_get_next(L,NULL); item!=NULL;) {
		tmp=item;
		item=list_get_next(L,item);
		if(tmp->p !=NULL)
			free(tmp->p);
		free(tmp);
	}
	free((void*)L);
}

void
list64_free_mem(list64 *L)
{
	list64_item *item,*tmp;
	if(L==NULL)
		return;
	for(item=list64_get_next(L,NULL); item!=NULL;) {
		tmp=item;
		item=list64_get_next(L,item);
		if(tmp->p !=NULL)
			free(tmp->p);
		free(tmp);
	}
	free((void*)L);
}



void
dump_list(list *L)
{
	list_item *tmp;

	printf("\n");
	for(tmp=L->l; tmp!=NULL; tmp=tmp->next) {
		printf("%d\t", tmp->val);
	}
	printf("\n");
}

void
dump_list64(list64 *L)
{
	list64_item *tmp;

	printf("\n");
	for(tmp=L->l; tmp!=NULL; tmp=tmp->next) {
		printf("%d\t", (int)tmp->val);
	}
	printf("\n");
}
