/************************************************************************
*
*  File name: role.c
*
*  Description: Role of nodes in a subgraph, related functions
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "list.h"
#include "hash.h"

/*************************** Global variables ****************************/

//Hash table to store all roles
Hash *Role_hash;

//translation array :
//from id and role to
role_st *Role_trans;

//Roles result table
Res_vec_tbl Roles_res_tbl;

//final results for role statistics
Role_res *Roles_final_res;

//Roles members 2 dim array [Role_id][node]
list **Roles_members;




//dictionary for roles id
//in Total: 30 roles id for motif size 30
//format: {motif id,role in the subgraph}
role_dic Role_dic[31]= {
	{0,0},
	{6,1},
	{6,2},
	{12,1},
	{12,2},
	{12,3},
	{14,1},
	{14,2},
	{14,3},
	{36,1},
	{36,2},
	{38,1},
	{38,2},
	{38,3},
	{46,1},
	{46,2},
	{74,1},
	{74,2},
	{74,3},
	{78,1},
	{78,2},
	{98,1},
	{102,1},
	{102,2},
	{102,3},
	{108,1},
	{108,2},
	{110,1},
	{110,2},
	{110,3},
	{238,1}
};





//An array of all possible 3-node ids and their roles
role_st role_3_arr[ROLE_3_ENTRIES] ={
	{ 6,{1,2,2} },
	{ 40,{2,1,2} },
	{ 192,{2,2,1} },

	{ 12,{2,1,3} },
	{ 34,{1,2,3} },
	{ 66,{2,3,1} },
	{ 96,{3,1,2} },
	{ 132,{1,3,2} },
	{ 136,{3,2,1} },

	{ 14,{1,2,3} },
	{ 42,{2,1,3} },
	{ 70,{1,3,2} },
	{ 168,{3,1,2} },
	{ 196,{2,3,1} },
	{ 224,{3,2,1} },

	{ 36,{1,1,2} },
	{ 72,{2,1,1} },
	{ 130,{1,2,1} },

	{ 38,{1,2,3} },
	{ 44,{2,1,3} },
	{ 104,{3,1,2} },
	{ 134,{1,3,2} },
	{ 194,{2,3,1} },
	{ 200,{3,2,1} },

	{ 46,{1,1,2} },
	{ 198,{1,2,1} },
	{ 232,{2,1,1} },

	{ 74,{2,3,1} },
	{ 76,{2,1,3} },
	{ 100,{3,1,2} },
	{ 138,{3,2,1} },
	{ 162,{1,2,3} },
	{ 164,{1,3,2} },

	{ 78,{1,2,2} },
	{ 170,{2,1,2} },
	{ 228,{2,2,1} },

	{ 98,{1,1,1} },
	{ 140,{1,1,1} },

	{ 102,{1,2,3} },
	{ 106,{3,1,2} },
	{ 142,{1,3,2} },
	{ 172,{2,1,3} },
	{ 204,{3,2,1} },
	{ 226,{2,3,1} },

	{ 108,{2,1,2} },
	{ 166,{1,2,2} },
	{ 202,{2,2,1} },

	{ 110,{2,1,3} },
	{ 174,{1,2,3} },
	{ 206,{2,3,1} },
	{ 230,{1,3,2} },
	{ 234,{3,2,1} },
	{ 236,{3,1,2} },

	{ 238,{1,1,1} }
};


/******************************* Externs *********************************/
extern Network *G_N;



/******************************* Functions *******************************/

//init translation array
//translate from
int
init_role_trans(int mtf_sz)
{
	int rc=RC_OK;
	int i;
	int mtf_id;

	//init dic
	Role_trans=(role_st*)calloc(256,sizeof(role_st));
	for(i=0;i<256;i++) {
		Role_trans[i].id=0;
		Role_trans[i].roles[0]=0;
		Role_trans[i].roles[1]=0;
		Role_trans[i].roles[2]=0;
	}

	if(mtf_sz == 3) {
		for(i=0;i<=30;i++){
			//motif id
			mtf_id=(int)Role_dic[i].id;
			Role_trans[mtf_id].id=mtf_id;
			//role gets the role ID
			Role_trans[mtf_id].roles[Role_dic[i].role]=(char)i;
		}
	}
	return rc;
}


void
free_role_trans()
{
	free(Role_trans);
}



int
init_role_hash(int mtf_sz)
{
	int rc=RC_OK;
	int i,j;
	role_st *role;

	//init hash
	Role_hash=(Hash*)calloc(1,sizeof(Hash));
	hash_init(Role_hash,HASH_SIZE,1);

	if(mtf_sz == 3) {
		for(i=0;i<ROLE_3_ENTRIES;i++){
			role=(role_st*)calloc(1,sizeof(role_st));
			role->id=role_3_arr[i].id;
			for(j=0;j<mtf_sz;j++)
				role->roles[j]=role_3_arr[i].roles[j];
			//insert to hash
			//nadav -added casting check..
			hash_insert(Role_hash,(void*)role,(int*)&role->id);
		}
	}
	return rc;
}




void
free_mem_role_hash()
{
	int i;
	list_item *item,*tmp;
	role_st *role;

	for(i=0;i<Role_hash->size;i++){
		if(Role_hash->h[i]!=NULL) {
			for(item=list_get_next(Role_hash->h[i],NULL); item!=NULL;) {
				tmp=item;
				item=list_get_next(Role_hash->h[i],item);
				if(tmp->p !=NULL){
					role=(role_st*)(tmp->p);
					free(role);
				}
				free(tmp);
			}
		}
	}
	free((void*)Role_hash->h);
}

int
init_roles_res_tbl(int rnd_net_num)
{
	int rc=RC_OK;
	int i;

	Roles_res_tbl.real=(int*)calloc(TOTAL_ROLES_3_NUM+1,sizeof(int));
	Roles_res_tbl.rand=(int**)calloc(rnd_net_num,sizeof(int*));
	for(i=0;i<rnd_net_num;i++)
		Roles_res_tbl.rand[i]=(int*)calloc(TOTAL_ROLES_3_NUM+1,sizeof(int));
	return rc;
}


int
get_roles_num_for_mtf_id(int64 id, int mtf_sz)
{
	if(mtf_sz!=3)
		return -1;
	switch((int)id){
	case 98:
	case 238:
		return 1;
	case 6:
	case 36:
	case 46:
	case 78:
	case 108:
		return 2;
	case 12:
	case 14:
	case 38:
	case 74:
	case 102:
	case 110:
		return 3;
	default:
		return -1;
		break;
	}
}




int
init_roles_members(int mtf_sz)
{
	int i,rc=RC_OK;
	if(mtf_sz==3) {
		//nadav - chnaged casting from int** to list **
		Roles_members=(list**)calloc(TOTAL_ROLES_3_NUM+1,sizeof(list*));
		for(i=1;i<=TOTAL_ROLES_3_NUM;i++)
			list_init(&Roles_members[i]);
	}else
		return RC_ERR;
	return rc;
}

void
free_role_members()
{
	int i;
	for(i=1;i<=TOTAL_ROLES_3_NUM;i++)
		list_free_mem(Roles_members[i]);
	free(Roles_members);
}




void
free_roles_res_tbl(int rnd_net_num)
{
  //int rc=RC_OK;
	int i;

	free(Roles_res_tbl.real);
	for(i=0;i<rnd_net_num;i++)
	  free(Roles_res_tbl.rand[i]);
	free(Roles_res_tbl.rand);
}

//update the Role vector result, in Roles_res_tbl
//according to N->roles_vec
int
update_roles_res(char **net_roles_vec,int vertices_num, int *tbl_roles_vec)
{
	int rc=RC_OK;
	int i,j;
	for(i=1;i<=vertices_num;i++)
		for(j=0;j<=TOTAL_ROLES_3_NUM;j++){
			tbl_roles_vec[j]+=net_roles_vec[i][j];
		}
		return rc;
}

//fill Roles_members according to results of roles counts in the network
//tbl_roles_vec is an array of int , for each tbl_roles_vec[ROLE_ID]="number of members"
int
fill_roles_members(int nodes_num,char **roles_vec,int *tbl_roles_vec, list**roles_members)
{
	int rc=RC_OK;
	int i,j;
	for(i=1;i<=TOTAL_ROLES_3_NUM;i++) {
		for(j=1;j<=nodes_num;j++){
			if(roles_vec[j][i]==1)
				list_insert(roles_members[i],j,NULL);
		}
		//check that the number of members for that role is equal to
		//the number of tbl_roles_vec
		if(roles_members[i]->size != tbl_roles_vec[i])
			printf("num of members of role %d is wrong\n",i);
	}
	return rc;
}


void
calc_roles_indeg_outdeg_stats(int *indeg,int *outdeg)
{
	int j;
	list_item *l_i;
	double indeg_sum=0,indeg_x;
	double indeg_max;
	double indeg_sigma=0;
	double indeg_power;
	double outdeg_sum=0,outdeg_x;
	double outdeg_max;
	double outdeg_sigma=0;
	double outdeg_power;
	int node;


	//for all roles id
	for(j=1;j<=TOTAL_ROLES_3_NUM;j++) {
		// calc degree mean val
		indeg_sum=0;outdeg_sum=0;
		indeg_max=0;outdeg_max=0;
		indeg_sigma=0;outdeg_sigma=0;
		for(l_i=list_get_next(Roles_members[j],NULL);
		l_i!=NULL;l_i=list_get_next(Roles_members[j],l_i)) {
			node=l_i->val;
			indeg_sum+=indeg[node];
			outdeg_sum+=outdeg[node];
			if(indeg[node] > indeg_max)
				indeg_max=indeg[node];
			if(outdeg[node] > outdeg_max)
				outdeg_max=outdeg[node];
		}
		if(Roles_members[j]->size > 0) {
			Roles_final_res[j].indeg_mean = (double)indeg_sum / (double)Roles_members[j]->size;
			Roles_final_res[j].outdeg_mean = (double)outdeg_sum / (double)Roles_members[j]->size;
		} else {
			Roles_final_res[j].indeg_mean=0;
			Roles_final_res[j].outdeg_mean=0;
		}
		Roles_final_res[j].indeg_max=indeg_max;
		Roles_final_res[j].outdeg_max=outdeg_max;

		if( Roles_members[j]->size > 1) {
			indeg_power=0;outdeg_power=0;
			indeg_sigma=0;outdeg_sigma=0;
			//calc in/out degree standard deviation
			for(l_i=list_get_next(Roles_members[j],NULL);
			l_i!=NULL;l_i=list_get_next(Roles_members[j],l_i)) {
				node=l_i->val;

				//indeg sigma
				indeg_x=indeg[node];
				if( (indeg_x-Roles_final_res[j].indeg_mean) == 0)
					indeg_power = 0;
				else
					indeg_power=pow((double)(indeg_x - Roles_final_res[j].indeg_mean), (double)2);
				indeg_sigma += indeg_power/(double)(Roles_members[j]->size-1);

				//outdeg sigma
				outdeg_x=outdeg[node];
				if( (outdeg_x-Roles_final_res[j].outdeg_mean) == 0)
					outdeg_power = 0;
				else
					outdeg_power=pow((double)(outdeg_x - Roles_final_res[j].outdeg_mean), (double)2);
				outdeg_sigma += outdeg_power/(double)(Roles_members[j]->size-1);
			}
		}
		Roles_final_res[j].indeg_std_dev=(double)sqrt(indeg_sigma);
		Roles_final_res[j].outdeg_std_dev=(double)sqrt(outdeg_sigma);

	}
}






void
calc_roles_deviation(int rnd_net_num)
{
	int i,j;
	double counts_sum=0,x;
	double sigma=0;
	double power;

	//for all roles id
	for(j=1;j<=TOTAL_ROLES_3_NUM;j++) {
		// calc counts mean val
		counts_sum=0;
		sigma=0;
		for(i=0; i<rnd_net_num; i++)
			counts_sum+=Roles_res_tbl.rand[i][j];
		if(rnd_net_num > 0)
			Roles_final_res[j].rand_mean = (double)counts_sum / (double)rnd_net_num;
		else
			Roles_final_res[j].rand_mean=0;

		if(rnd_net_num >1) {
			power=0;
			sigma=0;
			//calc counts standard deviation
			for(i=0; i<rnd_net_num; i++) {
				x=Roles_res_tbl.rand[i][j];
				if( (x-Roles_final_res[j].rand_mean) == 0)
					power = 0;
				else
					power=pow((double)(x - Roles_final_res[j].rand_mean), (double)2);
				sigma += power/(double)(rnd_net_num-1);
			}
		}
		Roles_final_res[j].rand_std_dev=(double)sqrt(sigma);
	}

}

void
calc_roles_pval(int rnd_net_num)
{
	int i,j;
	int real_pval=0;

	//for all roles id
	for(j=1;j<=TOTAL_ROLES_3_NUM;j++) {
		//calc real pval
		real_pval=0;
		for(i=0; i<rnd_net_num; i++) {
			if(Roles_final_res[j].real_count <= Roles_res_tbl.rand[i][j])
					real_pval++;
		}

		if(rnd_net_num>0)
			Roles_final_res[j].real_pval = (double)real_pval / (double)rnd_net_num;
		else
			Roles_final_res[j].real_pval=0;
	}
}

void
calc_roles_zscore()
{
	int j;

	//for all roles id
	for(j=1;j<=TOTAL_ROLES_3_NUM;j++) {
		if (Roles_final_res[j].rand_std_dev > 0) {
			Roles_final_res[j].real_zscore = (double)(Roles_final_res[j].real_count - Roles_final_res[j].rand_mean)
			/ (double)Roles_final_res[j].rand_std_dev;
		} else {
			Roles_final_res[j].real_zscore = 0;
		}
	}
}


void
calc_roles_ratio_deviation(int rnd_net_num)
{
	int i,j;
	double counts_sum=0,x;
	double sigma=0;
	double power;
	Motif_res *mtf_res=NULL;

	//for all roles id
	for(j=1;j<=TOTAL_ROLES_3_NUM;j++) {
		// calc ratio mean val
		counts_sum=0;
		for(i=0; i<rnd_net_num; i++) {
			counts_sum+=Roles_res_tbl.rand[i][j];
			mtf_res=(Motif_res*)Roles_final_res[j].mtf_res;
		}
		if((rnd_net_num > 0) && (mtf_res->rand_mean !=0)) {
			Roles_final_res[j].rand_ratio_mean = mtf_res->rand_mean/((double)counts_sum / (double)rnd_net_num);
			Roles_final_res[j].rand_ratio_mean *= Roles_final_res[j].role_app_inside_motif;
		}
		else
			Roles_final_res[j].rand_ratio_mean=0;

		if(rnd_net_num >1) {
			power=0;
			sigma=0;
			//calc counts standard deviation
			for(i=0; i<rnd_net_num; i++) {
				if(mtf_res->rand_mean!=0) {
					x=mtf_res->rand_mean/Roles_res_tbl.rand[i][j];
					x*=Roles_final_res[j].role_app_inside_motif;
					if( (x-Roles_final_res[j].rand_ratio_mean) == 0)
						power = 0;
					else
						power=pow((double)(x - Roles_final_res[j].rand_ratio_mean), (double)2);
				}else{
					power=0;
				}
				sigma += power/(double)(rnd_net_num-1);
			}
			Roles_final_res[j].rand_ratio_std_dev=(double)sqrt(sigma);
		}
	}
}

void
calc_roles_ratio_pval(int rnd_net_num)
{
	int i,j;
	int real_ratio_pval=0;
	Motif_res *mtf_res;
	double rand_ratio;

	//for all roles id
	for(j=1;j<=TOTAL_ROLES_3_NUM;j++) {
		//calc real pval
		real_ratio_pval=0;
		mtf_res=(Motif_res*)Roles_final_res[j].mtf_res;
		for(i=0; i<rnd_net_num; i++) {
			rand_ratio=mtf_res->rand_mean / Roles_res_tbl.rand[i][j];
			rand_ratio*=Roles_final_res[j].role_app_inside_motif;
			if(Roles_final_res[j].real_ratio <= rand_ratio)
					real_ratio_pval++;
		}

		if(rnd_net_num>0)
			Roles_final_res[j].real_ratio_pval = (double)real_ratio_pval / (double)rnd_net_num;
		else
			Roles_final_res[j].real_ratio_pval=0;
	}
}

void
calc_roles_ratio_zscore()
{
	int j;

	//for all roles id
	for(j=1;j<=TOTAL_ROLES_3_NUM;j++) {
		if (Roles_final_res[j].rand_ratio_std_dev > 0) {
			Roles_final_res[j].real_ratio_zscore = (double)(Roles_final_res[j].real_ratio - Roles_final_res[j].rand_ratio_mean)
			/ (double)Roles_final_res[j].rand_ratio_std_dev;
		} else {
			Roles_final_res[j].real_ratio_zscore = 0;
		}
	}
}



//calc Roles final results
int
calc_roles_final_res(int rnd_net_num, list64 *mtf_final_res)
{
	int rc=RC_OK;
	int i,j;
	Motif_res *mtf_res=NULL;
	list64_item *l_i;
	role_st *role_hash_entry;

	//allocated in search real network
	//Roles_final_res=(Role_res*)calloc(TOTAL_ROLES_3_NUM+1,sizeof(Role_res));
	for(j=1;j<=TOTAL_ROLES_3_NUM;j++) {
		Roles_final_res[j].mtf_id=Role_dic[j].id;
		//get pointer to Motif_res struct from final res
		if( (l_i=list64_get(mtf_final_res, (int64)Roles_final_res[j].mtf_id))==NULL)
			printf("Roles_final_res Error in getting mtf_res \n");
		else
			//update the pointer
			mtf_res=(Motif_res*)l_i->p;
			Roles_final_res[j].mtf_res=l_i->p;

		Roles_final_res[j].role_id=j;
		Roles_final_res[j].role=Role_dic[j].role;
		Roles_final_res[j].real_count=(double)Roles_res_tbl.real[j];
		//updsate role appearances in motif field
		for(i=0;i<3;i++) {
			//nadav - I added casting here, check if ok
			// and casting int64 to int
			role_hash_entry = (role_st*)hash_get(Role_hash,(int)mtf_res->id);
			if(Roles_final_res[j].role==role_hash_entry->roles[i])
				Roles_final_res[j].role_app_inside_motif++;
		}
		//real ratio : Motif Nreal / Role Nreal
		if(mtf_res->real_count > 0) {
			Roles_final_res[j].real_ratio=(double)mtf_res->real_count/(double)Roles_final_res[j].real_count;
			Roles_final_res[j].real_ratio*= Roles_final_res[j].role_app_inside_motif;
		}
		else
			Roles_final_res[j].real_ratio=0;

	}
	calc_roles_deviation(rnd_net_num);
	calc_roles_pval(rnd_net_num);
	calc_roles_zscore();

	calc_roles_ratio_deviation(rnd_net_num);
	calc_roles_ratio_pval(rnd_net_num);
	calc_roles_ratio_zscore();

	//calc indeg out deg statistics
	calc_roles_indeg_outdeg_stats(G_N->indeg,G_N->outdeg);
	return rc;
}







