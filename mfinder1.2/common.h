/************************************************************************
*
*  File name: common.h
*
*  Description: Header file
*
*  Copyright © 2002-2004 Weizmann Institute of Science,
*			   76100 Rehovot Israel, All rights reserved
*
*************************************************************************/
#ifndef __COMMON_H
#define __COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "bits.h"
#include "list.h"
#include "mat.h"
#include "role.h"
#include "random.h"


/******************* Definitions ************************/


#define VERSION         1.2

//general
#define TRUE        1
#define FALSE       0


// error codes
#define RC_ERR      -1
#define RC_OK        0

#define INFINITY     0x0fffffff
#define UNDEFINED    888888


#define REAL_NET     1
#define RAND_NET     2


#define SINGLE_EDGE   0
#define DOUBLE_EDGE   1



//general algorithm definitions
#define DEFAULT_MTF_SZ           3   //default motif size
#define DEFAULT_RAND_NETWORK_NUM 100 //default ranfdom entworks no.
#define ZSCORE_TH				2    //z-score threshold
#define PVAL_TH					0.01 // p-value threshold
#define MFACTOR_TH				1.1  // mfactor threshold
#define UNIQUE_TH				4    //uniqueness th
#define R_SWITCH_FACTOR      100


//Random network generating despair ratio
//this value used to calculate the maximum tries to find proper candidates to exchnage
//despair num=switchs_num*DESPAIR_RATIO
#define DESPAIR_RATIO     50
#define MAX_MEMBERS_LIST_SZ 1000   //maximal members list size
#define TOP_MOTIF_LIST_SZ   7   //default top motif list size
#define MAX_FULL_MOTIF_LIST_SZ 10
#define SRC_TRG_FORMAT 1
#define TRG_SRC_FORMAT 2

//rqandom network geberation params
#define ADD_EDGES_DISPAIR_RATIO     10
#define GEN_NETWORK_DISPAIR 10000

/******************* Structures *************************/



// Edge structure
typedef struct {
	int s;   //source vertex
	int t;   //target vertex
	int weight; //color val
}Edge;


//Network structure
typedef struct {
	char *name;
	int vertices_num;  //total number of vertices
	int edges_num;     //total number of edges

	Mat *mat;            //Matrix respresents the network
	Edge *e_arr;      //edges array
	int e_sin_num;    //single edges num
	Edge *e_arr_sin;  //single edges array
	int e_dbl_num;    //double edges num
	Edge *e_arr_dbl;  //double edges array
	int *indeg;  //in degree array (index i related node i)
	int *outdeg; //out degree array (index i related node i)
	int *doubledeg; //double edges degree array (index i related node i)
	int *self_edge; // entry i is true if node i has a self edge
	double *cluster_degree; // clustering coefficient series
	char **roles_vec;
	Rnd_pd *rnd_e_pd; //random probablity distridution struture ,for new prob app

	//Network stats
	int in_hub_deg; //ingoing hub deg
	int out_hub_deg; //outgoing hub deg
	int in_hub; //ingoing hub
	int out_hub; //outgoing
	int roots_num;
	int leafs_num;
	int con_vertices_num; //total connected vertices num

	//for probabilistic
	int hub; //hub node (node with max in+out deg)
	int hub_deg; //max in+out deg
	int *hub_edges_vec; //for prob alg efficiency
	int *hub_edges_indices; //for prob alg efficiency
	Mat *e_map; //mapping og all edges [i,j] to their edge index in e_arr

	//Layered network
	int num_of_layers; //num of layers in the network
	int layers_range; //range of edges (between layers)
	int *layer_node_num;//num of nodes in each layer starting at layer 1
	int *layer_edges_num;//num of edges from each layer starting from layer 1
}Network;




// Result table
typedef struct {
	list64 *real;
	list64 **rand_arr;
}Res_tbl;


//Motif info (Subgraph info) structure
typedef struct {
	int64 id;
	double count;
	double prob_count; //sigma (Ri/pi)
	double conc; //concentration in mili motifs (counts/total motif counts)*1000
	int hits; //hits num in prob approach
	list *members;  //for uniqueness
	list *all_members; //for list of all members of this motif
	int numberOfSelfEdges; // saves the number of self edges in the motif
	double conv_grade; //convergness grade
}Motif;

typedef struct {
	time_t start; //start time;
	time_t end; //end time
	double elapsed; //elapsed time
} Runtime;

typedef struct {
	int node;
	unsigned char role;
}Member;



//Statistics about random networks
typedef struct {
	double grass_reff; //effective number of randsom network in grassberger
	double clones_mean; //mean number of clones in grassberger alg
	double clones_stddev; //stddev of clones in grassberger alg
	double switch_ratio_mean; //ratio between success switches to total switches
	double switch_ratio_stddev;
}RndNetStat;



// The general program info and flags
typedef struct {
	int mtf_sz;   //motif size
	int rnd_net_num;  //number of random networks to generate

	int calc_unique_flag; //claculate unique appearances flag
	int undirected_flag; //input network is undirected graph
	int input_net_format; //format of inptu file SRC_TRG_FORMAT or TRG_SRC_FORMAT
	int calc_self_edges; //calculate motifs with self edges.
	int calc_weights;    //calculate motifs with weights.
	int list_members; //list members of all members of each motif
	int specific_subgraph_members;//list members of a specific subgraphs only
	int calc_roles; //calculate roles stats

	int use_metropolis;  //randomize networks with metropolis algorithm
	int use_stubs_method;  //randomize networks with ron's method
	int use_clustering;	   // randomize networks preserving clustering sequence

	double t_init;    // initial temperature for metropolis randomization
	double iteration_factor; // how many metropolis steps
	double e_thresh; // energy stop threshold for metropolis randomization


	int run_prob_app; //run probabilistic approach flag
	int efc_prob; //efficient prob algorithm
	int prob_base_samples_num; //probabilistic approach base samples number
	int *prob_total_samples_num; //probabilistic approach total samples number array
	int prob_converge_mode; //use convergness to decide the required no of samples
	double prob_conv_conc_th; //threshold concentration for convergness
	double prob_conv_diff; //avg diff requirement for convergness

	double mfactor_th; //mfactor to take into acount for th
	double zfactor_th; //zscore to take into acount for th
	int unique_th; //uniqueness threshold
	int force_unique_th; //force uniqueness threshold
	double pval_th; //pval threshold

	int quiet_mode; //quiet mode flag
	int dont_die; //dont die mode (wait to user interaction to finish)

	int dont_search_real;//do not search real network
	double r_switch_factor; //number of edges switch factor, when genrating random network
	int r_global_switch_mode; //Global switch :dont ignore unseccessfull switches in the
	int r_grassberger; //generate random network by grassberger alg
	int r_grass_colony_sz; //sizeof colonies in grassbrger alg
	int r_grass_colony_max_population; //ratio between colony size and max colony size
	int actual_n_eff; // Neff to use (1-5)
	int r_dont_conserve_mutuals; //do not conserve mutual edges

	int r_conserve_layers; //conserve layers in random netork
	int r_layers_num; //num of layers in input network (if working in conserve layers in rand nets)
	int *r_layer_sz; //layer size of input network (if working in conserve layers in rand nets)
	char out_fname[1024]; //output file name
	char log_fname[1024]; //log file name
	char mat_s_fname[1024]; //matlab short output file name
	char mat_l_fname[1024]; //matlab short output file name
	char mat_metrop_fname[1024]; //matlab metropolis output file name
	char inter_out_fname[1024]; //intermediate output file name
	char members_fname[1024]; //members of all motifs output file
	char roles_fname[1024]; //roles results output file name
	char rand_all_fname[1024]; //random network all results mat file
	char clust_fname[1024]; //clustering coefficient series file
	int max_members_list_sz; //maximum lengh of members list
	int top_motifs_list_sz; //top motifs list size
	int long_out_flag;  //short_output_flag
	int out_log_flag;   //output log file flag
	int out_s_mat_flag; //ouput matlab matrix short format
	int out_l_mat_flag; //ouput matlab matrix long format
	int out_s_c_mat_flag; //ouput matlab matrix short format in conc
	int out_metrop_mat_flag; //ouput metropolis format
	int out_intermediate; //ouput intermediate results
	int out_members; //output members list
	int out_roles; //output roles res list
	int out_rand_mat_flag; //output matlab matrix files of all subgraph counts
						   //for each random network
	int out_all_rand_networks; //output all random networks to text files
	int out_non_dangling_motifs; //ouput a list of all non-dangling motifs
	int out_clustering;	// output clustering series

	FILE *out_fp;
	FILE *log_fp;
	FILE *mat_s_fp;
	FILE *mat_l_fp;
	FILE *mat_metrop_fp;
	FILE *inter_out_fp;
	FILE *members_fp;
	FILE *roles_fp;
	FILE *mat_rand_fp;
	FILE *clust_fp;

	Runtime total_time;
	Runtime real_net_time;
	Runtime rand_net_time;
	RndNetStat rnstat;
}Gnrl_st;


// Statistics results about each subgraph
typedef struct {
	int64 id; //motif id
	int size;
	double real_count;  //real network motif count
	double real_pval;  //p value for real count
	double real_zscore;  //real count statistics distance
	int hits_num; //number of hits in probabilistic approach
	double conc_real; //real network motif concentration
	double conc_real_pval; //real network motif concentration pval
	double conc_real_zscore; //real network motif concentration zscore
	//Neff_st *eff_arr;
	int64 unique_appear; //unique appearances
	double rand_mean;  //mean count for random networks
	double rand_std_dev; //standard deviation for random networks
	double conc_rand_mean;  //mean concentraion of motif for random networks
	double conc_rand_std_dev; //standard deviation for motif concentration random networks
	list *all_members; //for list of all members of this motif
	int *all_rand_counts; //array of all rand counts
	double conv_grade; //convergness grade
	int dangling;//have dangling edges or not
} Motif_res;



/******************* Macros *****************************/




#define MAT(N,i,j)      ( *(char*)(N->mat->full.m+((i-1)*(N->vertices_num))+(j-1)) )



/******************* Prototypes *************************/


void
time_measure_start(Runtime *runtime);
int
count_subgraphs(Network *N, int mtf_sz, list64 **res_p, int net_type);
void
list_free_members_list(list *l);
void
list_free_all_members_list(list *l);
int
check_if_unique(Motif *mtf,int *mtf_vrtx_arr,int mtf_sz);
int
check_if_new_members(Motif *mtf,Member *mtf_vrtx_arr,int mtf_sz);
void
free_network_mem(Network *N);
void
free_subset_matrix_mem(Network *N);
int64
get_sn_motif_id(Network *SN);
int
motifs_search_real(Network *N);
int
duplicate_network(Network *SRC, Network **TRG_p, char *trg_name);
int
count_subgraphs_size_n(Network *N, int mtf_sz,list64 **res_p,int net_type,int rand_net_indx);
int
allocate_network(Network *SRC, Network **TRG_p, char *trg_name);
int
update_network(Network *RN,Network *N);
char*
my_itoa(int val, int base);
int
at_exit(int rc);
#endif
