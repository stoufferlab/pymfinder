/************************************************************************
*
*  File name: hash.h
*
*  Description: Header file
*  
*  Copyright © 2002-2004 Weizmann Institute of Science, 
*			   76100 Rehovot Israel, All rights reserved 
*
*************************************************************************/
#ifndef __HASH_H
#define __HASH_H

#include <stdio.h>


/******************* Definitions ************************/


//The following were chosen arbitrarily
#define HASH_SIZE      1013
#define HASH_KEY_BASE  17


/******************* Structures *************************/

typedef struct {
	int size;    // hash table size
	int bckt_sz; //backet array size
	list **h;
}Hash;


/******************* Prototypes *************************/

void
hash_init(Hash *hash,int size,int bckt_sz);

Hash*
multi_hash_init(int hash_tbl_num);

int
hash_get_key(Hash *hash,int *val);

void*
hash_get(Hash *hash, int val);

void
hash_insert(Hash *hash, void*item, int*val);

void
hash_free_mem(Hash *hash);

void
multi_hash_free_mem(Hash *multi_hash, int hash_tbl_num);





#endif



