/************************************************************************
*
*  File name: list.h
*
*  Description: Header file
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#ifndef __LIST_H
#define __LIST_H


#include <stdio.h>
#include "bits.h"



/******************* Definitions ************************/


/******************* Structures *************************/

// REGULAR list structure
//value field is 32 bit

typedef struct _list_item{
	int val;
	void *p;
	struct _list_item *next;
}list_item;

typedef struct{
	int size;
	list_item *l;
}list;

// EXTENDED list structure
//value field is 64 bit

typedef struct _list64_item{
	int64 val;
	void *p;
	struct _list64_item *next;
}list64_item;

typedef struct{
	int size;
	list64_item *l;
}list64;


/******************* Prototypes *************************/



// REGULAR list prototypes


void
list_init(list **L);

void
list_insert(list *L, int val, void *p);

int
list_delete(list *L, int val);

list_item*
list_get(list *L, int val);

list_item*
list_get_by_indx(list *L, int indx);

list_item*
list_get_next(list *L, list_item *item);

void
list_free_mem(list *L);

void
dump_list(list *L);





// EXTENDED list (64 bits) prototypes
void
list64_init(list64 **L);

void
list64_insert(list64 *L, int64 val, void *p);

int
list64_delete(list64 *L, int64 val);

list64_item*
list64_get(list64 *L, int64 val);

list64_item*
list64_get_by_indx(list64 *L, int indx);

list64_item*
list64_get_next(list64 *L, list64_item *item);

void
list64_free_mem(list64 *L);

void
dump_list64(list64 *L);



#endif




