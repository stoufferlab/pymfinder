/************************************************************************
*
*  File name: results.h
*
*  Description: Header file
*  
*  Copyright © 2002-2004 Weizmann Institute of Science, 
*			   76100 Rehovot Israel, All rights reserved 
*
*************************************************************************/
#ifndef __RESULTS_H
#define __RESULTS_H


/******************* Definitions ************************/


/******************* Structures *************************/


/******************* Prototypes *************************/


void
init_res_tbl(Res_tbl *res);
void 
res_tbl_mem_free_single(list64 *L);
void 
res_tbl_mem_free(Res_tbl *res);
void 
final_res_free(list64 *motif_res_lis);
void
calc_stat(double*arr, int sz, double* mean_p, double *stddev_p);
void
join_subgraphs_res(list64 **res_p, int mtf_sz, int rand_net_indx);
void
update_global_res_tbl(list64 *global_res,list64 *curr_res,int rand_net_indx);
int
calc_final_results(Res_tbl*res_tbl, list64 **final_res_p, list64 **final_res_all_p, int rnd_net_num);
int
calc_final_results_grassberger(Res_tbl*res_tbl, int sub_flag, list64 *res_sub_mtf, list64 **final_res_p, list64 **final_res_all_p, int rnd_net_num, double *weights_arr);
int
calc_generalized_motifs_final_results(Res_tbl *g_res_tbl,Res_tbl *g_cmplx_res_tbl, list64 **g_res_p, list64 **g_cmplx_res_p, int rnd_net_num);



#endif
