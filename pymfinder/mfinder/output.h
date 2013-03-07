/************************************************************************
*
*  File name: output.h
*
*  Description: Header file
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#ifndef __OUTPUT_H
#define __OUTPUT_H

#include <stdio.h>
#include "bits.h"


/******************* Definitions ************************/

/******************* Structures *************************/

/******************* Prototypes *************************/


void
dump_network(FILE *fp,Network *N);
void
dump_motif_matrix(FILE *fp, int64 mtf_id);
void
dump_motifs_res(list64 *res, char *network_name);
void
dump_members_lists(FILE *fp,list64 *res);
void
dump_final_res(FILE *fp,list64 *res, list64* full_res, int rnd_net_num);
void
dump_all_ids_N_data_matrix(FILE *fp, list64* full_res);
void
dump_all_ids_C_data_matrix(list64* full_res);

#if 0
void
dump_all_ids_data_matrix_full(list64* full_res);
#endif

void
dump_time_measure(FILE* fp,char *text, Runtime *runtime);
void
print_public_help();
void
print_private_help();
int
output_results(list64 *final_res, list64* final_res_all);
int
dump_roles_results(FILE *fp);
void
output_network_to_text_file(Network *N, char *fname);






#endif
