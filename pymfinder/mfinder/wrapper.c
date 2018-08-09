/************************************************************************
*
*
*
*
*
*
*
*
*************************************************************************/

#include "globals.h"
#include "mat.h"
#include "motif_ids.h"
#include "wrapper.h"

/*************************** Global variables ***************************/

/******************************* Externs ********************************/
// input file name
extern char *input_network_fname;

// crazy wicked wild "general program info and flags"
extern Gnrl_st GNRL_ST;

// confusing as hell result table
extern Res_tbl RES_TBL;

// a somewhat pared-down version of the results
extern list64 *final_res, *final_res_all;

// from the globals in main.c
extern int DEBUG_LEVEL;

// we need the network as an extern
extern Network *G_N;

/************************************************************************/

void set_default_options(){
  //copied from process_input_args in main.c
  GNRL_ST.mtf_sz=DEFAULT_MTF_SZ;
  GNRL_ST.rnd_net_num=DEFAULT_RAND_NETWORK_NUM;
  GNRL_ST.t_init=T0;
  GNRL_ST.iteration_factor=ITERATION_F;
  GNRL_ST.e_thresh=ETHRESH;
  GNRL_ST.use_stubs_method=FALSE;
  GNRL_ST.long_out_flag=FALSE;
  GNRL_ST.calc_unique_flag=TRUE;
  GNRL_ST.calc_roles=FALSE;
  GNRL_ST.input_net_format=SRC_TRG_FORMAT;
  GNRL_ST.undirected_flag=FALSE;
  GNRL_ST.calc_self_edges=FALSE;
  GNRL_ST.calc_weights=FALSE;
  GNRL_ST.run_prob_app=FALSE;
  GNRL_ST.prob_base_samples_num=0;
  GNRL_ST.prob_converge_mode=FALSE;
  GNRL_ST.prob_conv_diff=CONVERGNESS_DIFF_CONST;
  GNRL_ST.prob_conv_conc_th=CONC_THRESHOLD;
  GNRL_ST.unique_th=UNIQUE_TH;
  GNRL_ST.force_unique_th=FALSE;
  GNRL_ST.mfactor_th=MFACTOR_TH;
  GNRL_ST.zfactor_th=ZSCORE_TH;
  GNRL_ST.pval_th=PVAL_TH;
  GNRL_ST.max_members_list_sz=1000;
  GNRL_ST.top_motifs_list_sz=TOP_MOTIF_LIST_SZ;
  GNRL_ST.out_log_flag=FALSE;
  GNRL_ST.out_s_mat_flag=FALSE;
  GNRL_ST.out_l_mat_flag=FALSE;
  GNRL_ST.out_s_c_mat_flag=FALSE;
  GNRL_ST.out_members=FALSE;
  GNRL_ST.out_rand_mat_flag=FALSE;
  GNRL_ST.out_roles=FALSE;
  GNRL_ST.out_clustering=FALSE;
  GNRL_ST.quiet_mode=FALSE;
  GNRL_ST.use_metropolis=FALSE;
  GNRL_ST.use_clustering=FALSE;
  GNRL_ST.dont_search_real=FALSE;
  GNRL_ST.out_intermediate=FALSE;
  GNRL_ST.out_non_dangling_motifs=FALSE;
  GNRL_ST.r_switch_factor=(double)R_SWITCH_FACTOR;
  GNRL_ST.r_global_switch_mode=FALSE;
  GNRL_ST.list_members=FALSE;
  GNRL_ST.specific_subgraph_members=0; //no specific subgraph members
  GNRL_ST.efc_prob=FALSE;
  GNRL_ST.out_all_rand_networks=FALSE;
  GNRL_ST.r_grassberger=FALSE;
  GNRL_ST.r_grass_colony_sz=0;
  GNRL_ST.r_grass_colony_max_population=MAX_COLONY_SZ_RATIO;
  GNRL_ST.r_dont_conserve_mutuals=FALSE;
  GNRL_ST.dont_die=FALSE;
  GNRL_ST.r_conserve_layers=FALSE;
  GNRL_ST.r_layers_num=0;
}

int
load_network_from_array(Network **N_p, int* edges, int edges_num)
{
  int   rc=RC_OK,s,t,i,j,k,weight;
  list_item *l_e;
  int max_node=0;
  //FILE *fp;
  Network *N;
  int e_idx;
  int self_edge_exists=FALSE;
  int matVal = 0;

  N=(Network*)calloc(1,sizeof(Network));
  N->name="I came from an array, sucka!";
  N->edges_num = 0;
  N->e_self_num = 0;

  for(i=1;i<=edges_num;i++){
    s = edges[3*(i-1)+1];
    t = edges[3*(i-1)+2];
    if (s==t){
      N->e_self_num += 1;
      if (GNRL_ST.calc_self_edges == TRUE)
	N->edges_num += 1;
    }else
      N->edges_num += 1;
  }

  N->e_arr=(Edge*)calloc(N->edges_num+1,sizeof(Edge));
  N->e_arr_self=(Edge*)calloc(N->e_self_num+1,sizeof(Edge));
  j = k = 0;
  for(i=1;i<=edges_num;i++){
    //SRC->TRG format
    if(GNRL_ST.input_net_format==SRC_TRG_FORMAT){
      s = edges[3*(i-1)+1];
      t = edges[3*(i-1)+2];
      weight = edges[3*(i-1)+3];
    }
    //TRG->SRC format
    else{
      t = edges[3*(i-1)+1];
      s = edges[3*(i-1)+2];
      weight = edges[3*(i-1)+3];
    }

    if (s==t){
      N->e_arr_self[++j].s=s;
      N->e_arr_self[j].t=t;
      N->e_arr_self[j].weight=weight;
      
      if(GNRL_ST.calc_self_edges == TRUE){
	N->e_arr[++k].s=s;
	N->e_arr[k].t=t;
	N->e_arr[k].weight=weight;

	if (s>max_node) max_node=s;
	if (t>max_node) max_node=t;
      }
    }else{
      N->e_arr[++k].s=s;
      N->e_arr[k].t=t;
      N->e_arr[k].weight=weight;

      if (s>max_node) max_node=s;
      if (t>max_node) max_node=t;
    }
  }    

  ///////////////////////////////////////////////////////////////////
  // ALL OF THE ITEMS BELOW ARE COPIED FROM load_network IN GLOBALS.H
  ///////////////////////////////////////////////////////////////////

  N->vertices_num=max_node;

  rc |= MatInit(&N->mat,N->vertices_num,SPARSE);
  if(rc == RC_ERR) {
    printf("Error : Memory allocation failure\n");
    at_exit(-1);
  }

  //assign matrix entries according to edges
  for(i=1; i<=N->edges_num; i++) {
    if(GNRL_ST.calc_weights == TRUE) {
      MatAsgn(N->mat,N->e_arr[i].s,N->e_arr[i].t,N->e_arr[i].weight);
    }else{
      if(MatGet(N->mat,N->e_arr[i].s,N->e_arr[i].t) != 0) {
	printf("Error : Duplicate appearance of the edge (%d,%d) \n\tfound in the Network\n ",
	       N->e_arr[i].s,N->e_arr[i].t);
	at_exit(-1);
      }else{
	MatAsgn(N->mat,N->e_arr[i].s,N->e_arr[i].t,1);
      }
    }
  }

  //if efc_prob app then
  //create e_map matrix mapping of edges [s,t] to their index in e_arr
  if(GNRL_ST.efc_prob){
    rc |= MatInit(&N->e_map,N->vertices_num,SPARSE);
    if(rc == RC_ERR) {
      printf("Error : Memory allocation failure\n");
      at_exit(-1);
    }
    for(i=1; i<=N->edges_num; i++)

      MatAsgn(N->e_map,N->e_arr[i].s,N->e_arr[i].t,i);
  }

  //check there are no Self edges- if there are any then output them to screen and  stop
  if(GNRL_ST.calc_self_edges == FALSE){
    for(i=1; i<=N->vertices_num; i++) {
      if(MatGet(N->mat,i,i)==1) {
	fprintf(stdout,"Self edges exist in Input Network!!!\n");
	fprintf(stdout,"Self Edges : (%d,%d)\t",i,i);
	self_edge_exists=TRUE;
      }
    }
    if(self_edge_exists==TRUE)
      at_exit(-1);
  }

  //allocate and fill arrays of single edges and double edges
  N->e_arr_sin=(Edge*)calloc(N->edges_num+1,sizeof(Edge));
  N->e_arr_dbl=(Edge*)calloc((N->edges_num+1),sizeof(Edge));
  N->e_sin_num=0;
  N->e_dbl_num=0;
  N->roots_num=0;
  N->leafs_num=0;
  //allocate indeg and out deg arrays
  N->indeg=(int*)calloc(N->vertices_num+2,sizeof(int));
  N->outdeg=(int*)calloc(N->vertices_num+2,sizeof(int));
  N->doubledeg=(int*)calloc(N->vertices_num+2,sizeof(int));
  if(GNRL_ST.calc_self_edges == TRUE){
    N->self_edge=(int*)calloc(N->vertices_num+2,sizeof(int));
  }
  //actually matrix is sparse anyway now
  if(N->mat->type == SPARSE) {
    for(i=1;i<=N->vertices_num;i++) {
      for(j=1;j<=N->vertices_num;j++) {
	if( (matVal = MatGet(N->mat,i,j))) {
	  //if an edge and is not self edge
	  if(i != j){
	    //if the twin edge exists
	    if(MatGet(N->mat,j,i)){
	      //not inserted yet
	      if(j>i) {
		//double edge- this way always the twin pair has indexes 2x-1,2x
		N->e_arr_dbl[++N->e_dbl_num].s=i;
		N->e_arr_dbl[N->e_dbl_num].t=j;
		N->e_arr_dbl[N->e_dbl_num].weight=matVal;
		N->e_arr_dbl[++N->e_dbl_num].s=j;
		N->e_arr_dbl[N->e_dbl_num].t=i;
		N->e_arr_dbl[N->e_dbl_num].weight=matVal;
	      }
	    }
	    else {
	      //single edge
	      N->e_arr_sin[++N->e_sin_num].s=i;
	      N->e_arr_sin[N->e_sin_num].t=j;
	      N->e_arr_sin[N->e_sin_num].weight=matVal;
	    }
	  }
	  else //self edge
	    {
	      N->e_arr_sin[++N->e_sin_num].s=i;
	      N->e_arr_sin[N->e_sin_num].t=j;
	      N->e_arr_sin[N->e_sin_num].weight=matVal;
	    }
	}
      }
    }

    //fill in deg and out deg arrays
    for(i=0;i<=N->vertices_num;i++) {
      N->indeg[i]=0;
      N->outdeg[i]=0;
      if(GNRL_ST.calc_self_edges == TRUE)
	{
	  N->self_edge[i]=FALSE;
	}
    }
    for (i=1; i<=N->vertices_num;i++){
      if(N->mat->spr->m[i].to==NULL)
	N->outdeg[i]=0;
      else
	N->outdeg[i]=N->mat->spr->m[i].to->size;
      if(N->mat->spr->m[i].from==NULL)
	N->indeg[i]=0;
      else
	N->indeg[i]=N->mat->spr->m[i].from->size;
      if((N->mat->spr->m[i].self_edge == 1) && (GNRL_ST.calc_self_edges == TRUE))
	{
	  N->self_edge[i] = TRUE;
	}
    }
  }

  //statistics and global info about the network
  N->con_vertices_num=0;
  N->hub_deg=0;
  N->hub=0;
  N->in_hub_deg=0;
  N->in_hub=0;
  N->out_hub_deg=0;
  N->out_hub=0;
  //calc total num of connected vertices
  //and find hub preferenced
  for(i=1; i<=N->vertices_num;i++){
    if( (N->indeg[i]!=0) || (N->outdeg[i]!=0) )
      N->con_vertices_num++;
    if( ((N->indeg[i] + N->outdeg[i]) > N->hub_deg) ){
      N->hub_deg=N->indeg[i] + N->outdeg[i];
      N->hub=i;
    }
    if( N->indeg[i] > N->in_hub_deg){
      N->in_hub_deg=N->indeg[i];
      N->in_hub=i;
    }
    if( N->outdeg[i] > N->out_hub_deg){
      N->out_hub_deg=N->outdeg[i];
      N->out_hub=i;
    }
  }

  //if calc roles then init array or roles vector for all nodes
  if(GNRL_ST.calc_roles==TRUE && GNRL_ST.mtf_sz==3) {
    if( (N->roles_vec=(char**)calloc(N->vertices_num+1,sizeof(char*))) ==NULL)
      return RC_ERR;
    for(i=1;i<=N->vertices_num;i++)
      if( (N->roles_vec[i]=(char*)calloc(TOTAL_ROLES_3_NUM+1,sizeof(char))) ==NULL)
	return RC_ERR;
  }
  //if use clusterring in random network genereation
  if(GNRL_ST.use_clustering){
    N->cluster_degree=(double*)calloc(N->vertices_num+2,sizeof(double));
    for(i=0;i<=N->vertices_num;i++) {
      N->cluster_degree[i]=0.0;
    }
    //fill clustering series
    clustering_series(N,NULL);
  }
  //if efc prob approach then init hub_edges_vec and hub_edges_indices
  if(GNRL_ST.efc_prob){
    N->hub_edges_vec=(int*)calloc(N->edges_num+1,sizeof(int));
    N->hub_edges_indices=(int*)calloc(N->hub_deg+1,sizeof(int));

    //fill hub vector and hub indices
    i=0;
    for(l_e=list_get_next(N->e_map->spr->m[N->hub].to,NULL);l_e!=NULL;
	l_e=list_get_next(N->e_map->spr->m[N->hub].to,l_e)) {
      e_idx=*(int*)l_e->p;
      N->hub_edges_indices[++i]=e_idx;
      N->hub_edges_vec[e_idx]=i;
    }
    for(l_e=list_get_next(N->e_map->spr->m[N->hub].from,NULL);l_e!=NULL;
	l_e=list_get_next(N->e_map->spr->m[N->hub].from,l_e)) {
      e_idx=*(int*)l_e->p;
      N->hub_edges_indices[++i]=e_idx;
      N->hub_edges_vec[e_idx]=i;
    }
    if(i!=N->hub_deg)
      printf("Error in Hub degree\n");
  }
  //if conserve layers in random netowrks
  if(GNRL_ST.r_conserve_layers==TRUE){
    N->num_of_layers=GNRL_ST.r_layers_num;
    N->layer_node_num=(int*)calloc(N->num_of_layers+1,sizeof(int));
    N->layer_edges_num=(int*)calloc(N->num_of_layers+1,sizeof(int));
    //copy num of nodes in each layer
    for(i=1;i<=N->num_of_layers;i++)
      N->layer_node_num[i]=GNRL_ST.r_layer_sz[i];
    //run through e_arr and fill N->layer_edges_num
    //this is required for the randomizing process - in order to swtich edges
    //between nodes in the same layers only
    max_node=0;
    j=1;
    for(i=1;i<=N->num_of_layers;i++){
      max_node+=N->layer_node_num[i];
      while( (N->e_arr[j].s<=max_node) && (j<=N->edges_num) ){
	N->layer_edges_num[i]++;
	j++;
      }
    }
    //sanity check - that num of edges in  layers sum to total num of edges
    j=0;
    for(i=1;i<=N->num_of_layers;i++)
      j+=N->layer_edges_num[i];
    if(j!=N->edges_num){
      printf("\nERROR: in '-rcl' flag, check layers info\n");
      at_exit(-1);
    }
  }
  if(DEBUG_LEVEL>11)
    dump_network(stdout,N);
  *N_p=N;
  return rc;
}

/*
* function : read_network
*   read network from different input types
*
* arguments:
*   N_p - reference to network structure
*   mfinderi - mfinder wrapper input struct
*
* return values:
*   RC_OK - if not error occured
*   RC_ERR - if error occured
*/
int
read_network(Network **N_p, mfinder_input mfinderi)
{
  int rc=RC_OK;

  // load network from input file
  if(mfinderi.Filename != NULL){
    rc = load_network(&G_N,mfinderi.Filename);
    if (rc == RC_ERR) {
      printf("load network from file failed\n");
      return RC_ERR;
    }
  }
  // load network from array of edge info
  // this exists for streamlined use of the python module
  else{
    rc = load_network_from_array(&G_N,mfinderi.Edges,mfinderi.NumEdges);
    if (rc == RC_ERR) {
      printf("load network from array failed\n");
      return RC_ERR;
    }
  }

  return rc;
}

/********************************************************
 * function : single_connected_component
 *  Given a motif id and size,
 *  Return whether there is a single connected component
 * arguments:
 *   id
 *   mtf_siz
 * return values:
 *   TRUE/FALSE
 *********************************************************/

int single_connected_component(int64 id,int mtf_sz){
  int i,j,k;
  int added;

  //fprintf(stderr,"checking the scc for id=%lli\n",id);

  // build the adjacency matrix
  Matrix *M;
  M = init_matrix(mtf_sz);
  fill_mat_id(M,id);

  // collect members of the connected component here
  int cc[mtf_sz+1];
  int checked[mtf_sz+1];
  for(i=1;i<=mtf_sz;i++){
    cc[i] = FALSE;
    checked[i] = FALSE;
  }

  // add the first node to the connected component
  cc[1] = TRUE;
  
  // iteratively attempt to add more rows/columns to the cc
  added = TRUE;
  while(added == TRUE){
    added = FALSE;
    for(i=1;i<=mtf_sz;i++)
      if(cc[i] == TRUE && checked[i] == FALSE){
        checked[i] = TRUE;
        for(j=1;j<=mtf_sz;j++)
          if(MTRX(M,i,j)==1 || MTRX(M,j,i)==1)
            if(cc[j] == FALSE){
              cc[j] = TRUE;
              added = TRUE;
            }
      }
  }

  free_matrix(M);
  free(M);     

  for(i=1;i<=mtf_sz;++i)
    if(cc[i] == FALSE)
      return FALSE;

  return TRUE;
}

/********************************************************
 * function : list_all_motifs
 *  List all motifs for a given size
 *  --> copied largely from calc_final_results in results.c
 * arguments:
 *   AAA
 *   BBB
 *   CCC
 * return values:
 *   XXX
 *   YYY
 *********************************************************/

list64* list_motifs(int mtf_sz){
  int i,j,illegal_id=FALSE;
  list64_item *l_tmp;
  int64 id,rep_id,mask_bit,mask_r,mask_c;

  list64* iso_list;
  list64* id_list_uniq;
  list64* id_list_scc_false;
  list64_init(&id_list_uniq);
  list64_init(&id_list_scc_false);

  for(id=0;id<=(int64)(pow(2,mtf_sz*mtf_sz)-1);id++){
      illegal_id=FALSE;
      
      /*
      //Check if we've already found this id via the motif isomorphisms
      if(list64_get(id_list_all,(int64)id)!=NULL)
        illegal_id = TRUE;
      if(illegal_id)
        continue;
      */

      //Check if the motif contains self edges
      for(i=0;i<mtf_sz;i++) {
        //bits on the diagonal
        mask_bit=(int64)pow(2,i*mtf_sz+i);
        if(id & mask_bit) {
          illegal_id = TRUE;
          break;
        }
      }
      if(illegal_id)
        continue;

      //Check if there is no isolated vertex (a vertex with no edges at all)
      //this is done by checking that there is at list one edge
      //at each row i or column i of the matrix
      //(by bit manipulation)
      for(i=0;i<mtf_sz;i++) {
        mask_r=0;
        mask_c=0;
        for(j=0;j<mtf_sz;j++) {
          mask_r |= (int64)pow(2,i*mtf_sz+j);
          mask_c |= (int64)pow(2,i+j*mtf_sz);
        }
        if( ! ((id & mask_r) || (id & mask_c)) ) {
          illegal_id = TRUE;
          break;
        }
      }
      if(illegal_id == TRUE)
        continue;

      // calculate the ids for all motif isomorphisms
      iso_list=calc_mtf_id_iso(id,mtf_sz);

      //Check if the isomorphism is null
      l_tmp=list64_get_next(iso_list,NULL);
      if(l_tmp==NULL)
        illegal_id = TRUE;
      if(illegal_id == TRUE)
        continue;

      //Check if the isomorphism is already in the list
      rep_id=l_tmp->val;
      if(list64_get(id_list_uniq,(int64)rep_id)!=NULL)
        illegal_id = TRUE;
      if(illegal_id == TRUE)
        continue;

      //Check if the motif has a single connected component
      if(single_connected_component(rep_id,mtf_sz) == FALSE)
        //list64_insert(id_list_scc_false,(int64)rep_id,NULL);
        illegal_id = TRUE;
      if(illegal_id == TRUE)
        continue;
              
      list64_insert(id_list_uniq,(int64)rep_id,NULL);
      
      /*
      for(l_tmp=list64_get_next(iso_list,NULL); l_tmp !=NULL; l_tmp=list64_get_next(iso_list,l_tmp)) {
        rrep_id=l_tmp->val;
        if(list64_get(id_list_all,(int64)rrep_id)==NULL){
          list64_insert(id_list_all,(int64)rrep_id,NULL);
          //fprintf(stdout,"iso_id %lli == iso_id %lli\n",rep_id,rrep_id);
        }
      }
      */

      list64_free_mem(iso_list);
  }

  list64_free_mem(id_list_scc_false);

  return id_list_uniq;
  //list64_free_mem(id_list_uniq);
}

list* motif_edges(int64 id,int mtf_sz){
  int i,j;
  Matrix *M;
  Edge *e;
  list *motif_edges;
  list_init(&motif_edges);

  M = init_matrix(mtf_sz);
  fill_mat_id(M,id);

  for(i=1;i<=mtf_sz;i++)
    for(j=1;j<=mtf_sz;j++)
      if(MTRX(M,i,j)==1){
        e=(Edge*)calloc(1,sizeof(Edge));
        e->s = i;
        e->t = j;
        e->weight = 1;

        list_insert(motif_edges,0,(void*)e);
      }

  free_matrix(M);
  free(M);

  return motif_edges;
}

/********************************************************
 * function : random_network
 *	Randomize an input network
 * arguments:
 *   AAA
 *   BBB
 *   CCC
 * return values:
 *   XXX
 *   YYY
 *********************************************************/

list* random_network(mfinder_input mfinderi){
  set_default_options();

  // turn on quiet mode
  GNRL_ST.quiet_mode=TRUE;

  // ignore self edges
  GNRL_ST.calc_self_edges=FALSE;

  // general initialization
  int rc = gnrl_init();
  if (rc == RC_ERR) {
    printf("general init failed\n");
    return NULL;
  }

  // initialize the random seed
  init_random_seed();

  // read in the network
  rc = read_network(&G_N,mfinderi);
  if(rc==RC_ERR){
    fprintf(stderr,"load network failed\n");
    return NULL;
  }

  // generate a new (to be randomized) network
  Network *N;

  // typical single/double switching algorithm
  if(mfinderi.UseMetropolis == 0){ // || GNRL_ST.mtf_sz <= 3){
    //fprintf(stderr,"randomizing in typical single-double fashion\n");
    double switch_ratio;
    rc = gen_rand_network_switches_method_conserve_double_edges(&N,&switch_ratio);
    if (rc == RC_ERR) {
      fprintf(stderr,"single-double randomize network failed\n");
      return NULL;
    }
    rc = update_network(N,G_N);
    if (rc == RC_ERR) {
      fprintf(stderr,"randomize network failed to maintain node statistics\n");
      return NULL;
    }
  }
  // metropolis algorithm to have same triad consensus
  else{
    //printf("randomize while using the metropolis algorithm to preserve triads\n");
    int j,real_vec13[14];
    Res_tbl met_res_tbl;

    for(j=1;j<=13;j++)
      real_vec13[j]=0;
    init_res_tbl(&met_res_tbl);

    met_motifs_search_real(G_N,&met_res_tbl,real_vec13);
    list64_free_mem(met_res_tbl.real);
    //res_tbl_mem_free(met_res_tbl);

    rc = gen_rand_network_metrop(&N,real_vec13);
    if (rc == RC_ERR) {
      fprintf(stderr,"metropolis randomize network failed\n");
      return NULL;
    }
    rc = update_network(N,G_N);
    if (rc == RC_ERR) {
      fprintf(stderr,"randomize network failed to maintain node statistics\n");
      return NULL;
    }
  }

  list *randomized_edges;
  list_init(&randomized_edges);
  
  int i;
  Edge *e;
  // add the self edges that were never randomized
  for(i=1;i<=N->e_self_num;++i){
    e=(Edge*)calloc(1,sizeof(Edge));
    e->s = N->e_arr_self[i].s;
    e->t = N->e_arr_self[i].t;
    e->weight = N->e_arr_self[i].weight;

    list_insert(randomized_edges,0,(void*)e);
  }
  // add the now randomized edges
  for(i=1;i<=N->edges_num;++i){
    e=(Edge*)calloc(1,sizeof(Edge));
    e->s = N->e_arr[i].s;
    e->t = N->e_arr[i].t;
    e->weight = N->e_arr[i].weight;

    list_insert(randomized_edges,0,(void*)e);
  }

  // release memory of various objects
  free_network_mem(N);
  free(N);

  return randomized_edges;
}

/********************************************************
 * function : motif_structure
 *	Calculate the motif structure statistics
 * arguments:
 *   AAA
 *   BBB
 *   CCC
 * return values:
 *   XXX
 *   YYY
 *********************************************************/

list64* motif_structure(mfinder_input mfinderi){
  set_default_options();

  // turn on quiet mode
  GNRL_ST.quiet_mode=TRUE;

  // what size motif are we talking about?
  GNRL_ST.mtf_sz = mfinderi.MotifSize;

  // against how many randomizations should we compare?
  GNRL_ST.rnd_net_num = mfinderi.NRandomizations;

  // ignore the uniqueness bit
  GNRL_ST.unique_th=0;
  GNRL_ST.calc_unique_flag=FALSE;

  // ignore self edges
  GNRL_ST.calc_self_edges=FALSE;

  // randomize with the metropolis algorithm to preserve the same triad consensus
  if(mfinderi.UseMetropolis == 1)
    GNRL_ST.use_metropolis=TRUE;

  // general initialization
  int rc = gnrl_init();
  if (rc == RC_ERR) {
    printf("general init failed\n");
    return NULL;
  }/*else
     printf("general init succeeded\n");*/

  // initialize the random seed
  init_random_seed();

  // read in the network
  rc = read_network(&G_N,mfinderi);
  if(rc==RC_ERR){
    printf("load network failed\n");
    return NULL;
  }

  // generate a fake network to search in
  Network *N;
  rc = duplicate_network(G_N,&N,"real_network");
  if (rc == RC_ERR) {
    printf("duplicate network failed\n");
    return NULL;
  }

  // exhaustive search motif size n
  rc = count_subgraphs(N, GNRL_ST.mtf_sz, &RES_TBL.real, REAL_NET);
  if (rc == RC_ERR) {
    printf("real search failed\n");
    return NULL;
  }

  // calc result after isomorphism of ids
  join_subgraphs_res(&RES_TBL.real, GNRL_ST.mtf_sz, 0);

  // conduct randomized network analysis
  if(GNRL_ST.rnd_net_num > 0){
    rc = process_rand_networks(&RES_TBL, GNRL_ST.mtf_sz);
    if (rc == RC_ERR) {
      printf("random network analysis failed\n");
      return NULL;
    }
  }

  //calculate final results
  calc_final_results(&RES_TBL, &final_res, &final_res_all, GNRL_ST.rnd_net_num);

  // release memory of various network objects
  free_network_mem(N);
  free(N);

  // release the results table memory
  res_tbl_mem_free(&RES_TBL);

  // return the crazy results table
  // per mfinder, this depends on the size of the motifs
  if(GNRL_ST.mtf_sz<=4){
    final_res_free(final_res);
    return final_res_all;
  }else{
    final_res_free(final_res);
    return final_res_all;
  }
}


/********************************************************
 * function : motif_participation
 *	Calculate the motif participation statistics
 * arguments:
 *   AAA
 *   BBB
 *   CCC
 * return values:
 *   XXX
 *   YYY
 *********************************************************/

list64* motif_participation(mfinder_input mfinderi){
  set_default_options();

  // turn on quiet mode
  GNRL_ST.quiet_mode=TRUE;

  // what size motif are we talking about?
  GNRL_ST.mtf_sz = mfinderi.MotifSize;

  // ignore the uniqueness bit
  GNRL_ST.unique_th=0;
  GNRL_ST.calc_unique_flag=FALSE;

  // ignore self edges
  GNRL_ST.calc_self_edges=FALSE;

  // keep track of motif members
  GNRL_ST.list_members=TRUE;
  GNRL_ST.max_members_list_sz=mfinderi.MaxMembersListSz;

  // general initialization
  int rc = gnrl_init();
  if (rc == RC_ERR) {
    printf("general init failed\n");
    return NULL;
  }

  // initialize the random seed
  init_random_seed();

  // read in the network
  rc = read_network(&G_N,mfinderi);
  if(rc==RC_ERR){
    printf("load network failed\n");
    return NULL;
  }

  // generate a fake real network or randomize the network before counting motifs and participation
  Network *N;
  // use the fake real network
  if (mfinderi.Randomize == 0){;
    rc = duplicate_network(G_N,&N,"real_network");
    if (rc == RC_ERR) {
      printf("duplicate network failed\n");
      return NULL;
    }
  }
  // randomize the network
  else{
    // typical single/double switching algorithm
    if(mfinderi.UseMetropolis == 0){ // || GNRL_ST.mtf_sz <= 3){
      //printf("randomizing in typical single-double fashion\n");
      double switch_ratio;
      rc = gen_rand_network_switches_method_conserve_double_edges(&N,&switch_ratio);
      if (rc == RC_ERR) {
	printf("single-double randomize network failed\n");
	return NULL;
      }
      rc = update_network(N,G_N);
      if (rc == RC_ERR) {
	printf("randomize network failed to maintain node statistics\n");
	return NULL;
      }
    }
    // metropolis algorithm to have same triad consensus
    else{
      //printf("randomize while using the metropolis algorithm to preserve triads\n");
      int j,real_vec13[14];
      Res_tbl met_res_tbl;

      for(j=1;j<=13;j++)
	real_vec13[j]=0;
      init_res_tbl(&met_res_tbl);

      met_motifs_search_real(G_N,&met_res_tbl,real_vec13);
      list64_free_mem(met_res_tbl.real);

      rc = gen_rand_network_metrop(&N,real_vec13);
      if (rc == RC_ERR) {
	printf("metropolis randomize network failed\n");
	return NULL;
      }
      rc = update_network(N,G_N);
      if (rc == RC_ERR) {
	printf("randomize network failed to maintain node statistics\n");
	return NULL;
      }
    }
  }

  //FILE *fp;
  //fp = fopen("fargus.net","w");
  //dump_network(stdout,N);
  //return NULL;
  //fclose(fp);

  //exhaustive search motif size n
  rc = count_subgraphs(N, GNRL_ST.mtf_sz, &RES_TBL.real, REAL_NET);
  if (rc == RC_ERR) {
    printf("real search failed\n");
    return NULL;
  }

  // release memory of various objects
  free_network_mem(N);
  free(N);

  //calc result after isomorphism of ids
  join_subgraphs_res(&RES_TBL.real, GNRL_ST.mtf_sz, 0);

  return RES_TBL.real;
}
