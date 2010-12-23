/************************************************************************
*
*  File name: roles.h
*
*  Description: Header file
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#ifndef __ROLE_H
#define __ROLE_H

#include "bits.h"
/******************* Definitions ************************/



//num of entries in role arr
//(equal to all possible ids of motifs 3 with symmetries
#define ROLE_3_ENTRIES       54
#define TOTAL_ROLES_3_NUM    30

/******************* Structures *************************/


typedef struct{
	int64 id;
	unsigned char roles[4];
}role_st;


typedef struct {
	int64 id;
	unsigned char role;
}role_dic;


typedef struct {
	int *real;
	int **rand;
}Res_vec_tbl;

typedef struct {
	int role_id; //role_id
	int64 mtf_id; //motif id
	void *mtf_res; //pointer to motif results (Motif_res*)
	unsigned char role; //role in id (1,2,3..)
	int role_app_inside_motif;  //number of appearances inside the subgraph (1 2 or 3)
	double real_count;  //real network Role count
	double real_pval;  //p value for Role count
	double real_zscore;  //real count statistics distance
	double rand_mean;  //mean count for random networks
	double rand_std_dev; //standard deviation for random networks

	double real_ratio;  //real network Role count
	double real_ratio_pval;  //p value for Role count
	double real_ratio_zscore;  //real count statistics distance
	double rand_ratio_mean;  //mean count for random networks
	double rand_ratio_std_dev; //standard deviation for random networks

	double indeg_mean; //in degree mean
	double indeg_max; //in degree max value
	double indeg_std_dev; //in degree standard deviation
	double outdeg_mean; //out degree mean
	double outdeg_max; //out degree max value
	double outdeg_std_dev; //out degree standard deviation
}Role_res;




/******************* Prototypes *************************/


//init Role_arr entries
//currently support 3 only
int init_role_hash(int mtf_sz);
void free_mem_role_hash();

//update role vec table accoding to roles vec of the network (real or rand)
int update_roles_res(char **net_roles_vec,int vertices_num,int *tbl_roles_vec);

int
calc_roles_final_res(int rnd_net_num, list64 *mtf_final_res);

//init role res table
//stores role results for real and rand networks
int init_roles_res_tbl();
void free_roles_res_tbl();
int fill_roles_members(int nodes_num,char **roles_vec,int *tbl_roles_vec, list** roles_members);
int init_roles_members(int mtf_sz);
void free_role_members();
int
init_role_trans(int mtf_sz);
//return number of roles for motif id
int
get_roles_num_for_mtf_id(int64 id, int mtf_sz);

#endif

