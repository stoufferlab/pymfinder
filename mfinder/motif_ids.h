/************************************************************************
*
*  File name: motif_ids.h
*
*  Description: Header file
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#ifndef __MOTIF_IDS_H
#define __MOTIF_IDS_H

#include <stdio.h>
#include "bits.h"
#include "list.h"
#include "mat.h"



/******************* Definitions ************************/



/******************* Structures *************************/

/******************* Prototypes *************************/
list64*
calc_mtf_id_iso(int64 id, int mtf_sz);
list64*
get_mtf_subid(int64 id, int mtf_sz);
void
fill_mat_id(Matrix *M, int64 id);
int64
get_mtf_subid_min_con_vrtx_less(int64 id, int mtf_sz);
int64
get_rep_mtf_id_3(int64 mtf_id);
void
fill_mat_id(Matrix *M, int64 id);
int
gen_subset_matrix(Network *N, Network **SN, list *vrtx_set, int mtf_sz);
int *
get_perm(int size,int index);

#endif
