/************************************************************************
 *
 * copied out of globals.h because I've changed the original functions
 *
 *
 *
 *
 *
 *
 *************************************************************************/
#include "common.h"
#include "stouffer.h"

/*************************** Global variables ****************************/


/******************************* Externs *********************************/

extern int DEBUG_LEVEL;
extern Gnrl_st GNRL_ST;

/******************************* Functions *******************************/

// free Network structure related allocated memory
void
free_network_mem(Network *N)
{
  int i;

  if(N!=NULL) {
    if(N->mat !=NULL){
      MatFree(N->mat);
      free(N->mat);
    }
    if(N->e_arr!=NULL)
      free(N->e_arr);
    if(N->e_arr_sin!=NULL)
      free(N->e_arr_sin);
    if(N->e_arr_dbl!=NULL)
      free(N->e_arr_dbl);
    if(N->e_arr_self!=NULL)
      free(N->e_arr_self);
    if(N->indeg!=NULL)
      free(N->indeg);
    if(N->outdeg!=NULL)
      free(N->outdeg);
    if(N->doubledeg!=NULL)
      free(N->doubledeg);
    if(N->cluster_degree!=NULL)
      free(N->cluster_degree);
    if(N->roles_vec!=NULL) {
      for(i=1;i<N->vertices_num;i++)
	free(N->roles_vec[i]);
      free(N->roles_vec);
    }
    if(GNRL_ST.efc_prob) {
      if(N->e_map!=NULL){
	MatFree(N->e_map);
	free(N->e_map);
      }
      if(N->hub_edges_indices!=NULL)
	free(N->hub_edges_indices);
      if(N->hub_edges_vec!=NULL)
	free(N->hub_edges_vec);
    }
    if(GNRL_ST.r_conserve_layers){
      free(N->layer_edges_num);
      free(N->layer_node_num);
    }
  }
}

/*
* function : load_network
*	Load network form input file
* arguments:
*   N_p - reference to network structure
*   network_fname - network file name
*
* return values:
*	RC_OK - if not error occured
*   RC_ERR - if error occured
*/
int
load_network(Network **N_p, char* network_fname)
{
  int   rc=RC_OK,s,t,i,j,weight;
  list *e_list, *e_self_list;
  list_item *l_e;
  int max_node=0;
  Edge *e;
  FILE *fp;
  Network *N;
  int e_idx;
  int self_edge_exists=FALSE;
  int matVal = 0;

  fp = fopen(network_fname,"rt");
  if(fp==NULL) {
    printf("\nError: Cannot open input file : %s\n\tF input file name and path\n",network_fname);
    at_exit(-1);
  }
  N=(Network*)calloc(1,sizeof(Network));
  N->name=network_fname;
  N->edges_num = 0;
  N->e_self_num = 0;

  //build edge array
  list_init(&e_list);
  list_init(&e_self_list);

  //SRC->TRG format
  if(GNRL_ST.input_net_format==SRC_TRG_FORMAT) {
    if(GNRL_ST.quiet_mode==FALSE)
      fprintf(stdout,"\tReading Network file in <Source,Target,Weight> Format\n");
    //read all edges one by one and insert each edge to edge list
    while ( (rc = fscanf( fp,"%d %d %d\n", &s, &t, &weight) != EOF) ) {
      //insert edge to edge list
      e=(Edge*)calloc(1,sizeof(Edge));
      e->s=s;
      e->t=t;
      if(GNRL_ST.calc_weights == TRUE)
	{
	  e->weight=weight;
	}
      else
	{
	  e->weight=1;
	}

      if(s == t)
	{
	  N->e_self_num++;
	  list_insert(e_self_list,0,(void*)e);
	  if(GNRL_ST.calc_self_edges == TRUE){
	    list_insert(e_list,0,(void*)e);
	    if (s>max_node) max_node=s;
	    if (t>max_node) max_node=t;
	    N->edges_num++;
	  }
	}
      else
	{
	  list_insert(e_list,0,(void*)e);
	  if (s>max_node) max_node=s;
	  if (t>max_node) max_node=t;
	  N->edges_num++;
	}
    }
  } else {
    //TRG->SRC fvormat
    if(GNRL_ST.quiet_mode==FALSE)
      fprintf(stdout,"\tReading Network file in <Target,Source,Weight> Format\n");
    //read all edges one by one and insert each edge to edge list
    while ( (rc = fscanf( fp,"%d %d %d\n", &t, &s, &weight) != EOF) ) {
      //insert edge to edge list
      e=(Edge*)calloc(1,sizeof(Edge));
      e->s=s;
      e->t=t;
      if(GNRL_ST.calc_weights == TRUE)
	{
	  e->weight=weight;
	}
      else
	{
	  e->weight=1;
	}

      if(s == t)
	{
	  N->e_self_num++;
	  list_insert(e_self_list,0,(void*)e);
	  if(GNRL_ST.calc_self_edges == TRUE){
	    list_insert(e_list,0,(void*)e);
	    if (s>max_node) max_node=s;
	    if (t>max_node) max_node=t;
	    N->edges_num++;
	  }
	}
      else
	{
	  list_insert(e_list,0,(void*)e);
	  if (s>max_node) max_node=s;
	  if (t>max_node) max_node=t;
	  N->edges_num++;
	}
    }
  }
  fclose(fp);

  //directed graph
  if(GNRL_ST.undirected_flag==FALSE) {
    //allocate e_arr fill it according to e_list and free e_list
    N->e_arr=(Edge*)calloc(N->edges_num+1,sizeof(Edge));

    for(i=N->edges_num,l_e=list_get_next(e_list,NULL);i>0;i--,l_e=list_get_next(e_list,l_e)){
      N->e_arr[i].s=((Edge*)(l_e->p))->s;
      N->e_arr[i].t=((Edge*)(l_e->p))->t;
      N->e_arr[i].weight=((Edge*)(l_e->p))->weight;
    }
  }
  //undirected graph
  else {
    N->edges_num = N->edges_num*2 - N->e_self_num;
    N->e_arr=(Edge*)calloc((N->edges_num+1),sizeof(Edge));
		
    for(i=N->edges_num,l_e=list_get_next(e_list,NULL);i>0;i--,l_e=list_get_next(e_list,l_e))
      {
	N->e_arr[i].s=((Edge*)(l_e->p))->s;
	N->e_arr[i].t=((Edge*)(l_e->p))->t;
	N->e_arr[i].weight=((Edge*)(l_e->p))->weight;
	if(((Edge*)(l_e->p))->s != ((Edge*)(l_e->p))->t) // not a self edge
	  {
	    i--;
	    N->e_arr[i].s=((Edge*)(l_e->p))->t;
	    N->e_arr[i].t=((Edge*)(l_e->p))->s;
	    N->e_arr[i].weight=((Edge*)(l_e->p))->weight;
	  }
      }
  }
  list_free_mem(e_list);

  //self edges
  //allocate e_arr_self and fill it according to e_self_list and free e_self_list
  N->e_arr_self=(Edge*)calloc(N->e_self_num+1,sizeof(Edge));
  for(i=N->e_self_num,l_e=list_get_next(e_self_list,NULL);i>0;i--,l_e=list_get_next(e_self_list,l_e)){
      N->e_arr_self[i].s=((Edge*)(l_e->p))->s;
      N->e_arr_self[i].t=((Edge*)(l_e->p))->t;
      N->e_arr_self[i].weight=((Edge*)(l_e->p))->weight;
    }
  list_free_mem(e_self_list);
 
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

/********************************************************
 * function : duplicate_network
 *	Duplicate
 * arguments:
 *   SRC - network source
 *   TRG - netowrk target
 *   network_fname - target network name
 *
 * return values:
 *	RC_OK - if not error occured
 *   RC_ERR - if error occured
 *********************************************************/
int
duplicate_network(Network *SRC, Network **TRG_p, char *trg_name)
{
  int i,j,rc=RC_OK;
  Network *TRG;
  list_item *l_e;
  TRG=(Network*)calloc(1,sizeof(Network));
  TRG->name=trg_name;
  TRG->vertices_num = SRC->vertices_num;
  TRG->edges_num = SRC->edges_num;
  TRG->e_self_num = SRC->e_self_num;
  MatInit(&TRG->mat, TRG->vertices_num,SPARSE);
  //FULL: TRG->m = (char*)malloc(TRG->vertices_num*TRG->vertices_num*sizeof(char));
  //init matrix according to input network
  if(SRC->mat->type==FULL) {
    for (i=1; i<=TRG->vertices_num; i++)
      for (j=1; j<=TRG->vertices_num; j++)
	MatAsgn(TRG->mat,i,j,MatGet(SRC->mat,i,j));
  } else {
    //sparse
    for (i=1; i<=TRG->vertices_num; i++) {
      for(l_e=list_get_next(SRC->mat->spr->m[i].to,NULL);l_e!=NULL;
	  l_e=list_get_next(SRC->mat->spr->m[i].to,l_e))
	MatAsgn(TRG->mat,i,l_e->val,*(int*)l_e->p);
      TRG->mat->spr->m[i].self_edge = SRC->mat->spr->m[i].self_edge;
    }
  }

  //copy edge arr from source network
  TRG->e_arr=(Edge*)calloc(TRG->edges_num+1,sizeof(Edge));
	
  for(i=1;i<=TRG->edges_num;i++) {
    TRG->e_arr[i].s=SRC->e_arr[i].s;
    TRG->e_arr[i].t=SRC->e_arr[i].t;
    TRG->e_arr[i].weight=SRC->e_arr[i].weight;
  }
  //copy e_arr_sin and e_arr_dup
  TRG->e_sin_num=SRC->e_sin_num;
  TRG->e_dbl_num=SRC->e_dbl_num;
  TRG->e_arr_sin=(Edge*)calloc((unsigned int)TRG->e_sin_num+1,sizeof(Edge));
  TRG->e_arr_dbl=(Edge*)calloc((unsigned int)TRG->e_dbl_num+1,sizeof(Edge));

  for(i=1;i<=TRG->e_sin_num;i++) {
    TRG->e_arr_sin[i].s=SRC->e_arr_sin[i].s;
    TRG->e_arr_sin[i].t=SRC->e_arr_sin[i].t;
    TRG->e_arr_sin[i].weight=SRC->e_arr_sin[i].weight;
  }
  for(i=1;i<=TRG->e_dbl_num;i++) {
    TRG->e_arr_dbl[i].s=SRC->e_arr_dbl[i].s;
    TRG->e_arr_dbl[i].t=SRC->e_arr_dbl[i].t;
    TRG->e_arr_dbl[i].weight=SRC->e_arr_dbl[i].weight;
  }

  //copy e_arr_self
  TRG->e_arr_self=(Edge*)calloc(TRG->e_self_num+1,sizeof(Edge));
  for(i=1;i<=TRG->e_self_num;i++) {
      TRG->e_arr_self[i].s=SRC->e_arr_self[i].s;
      TRG->e_arr_self[i].t=SRC->e_arr_self[i].t;
      TRG->e_arr_self[i].weight=SRC->e_arr_self[i].weight;
    }

  TRG->indeg=(int*)calloc((unsigned int)TRG->vertices_num+2,sizeof(int));
  TRG->outdeg=(int*)calloc((unsigned int)TRG->vertices_num+2,sizeof(int));
  TRG->doubledeg=(int*)calloc((unsigned int)TRG->vertices_num+2,sizeof(int));
  for(i=1;i<=TRG->vertices_num+1;i++) {
    TRG->indeg[i]=SRC->indeg[i];
    TRG->outdeg[i]=SRC->outdeg[i];
    TRG->doubledeg[i]=SRC->doubledeg[i];
  }
  if(GNRL_ST.use_clustering){
    TRG->cluster_degree=(double*)calloc(TRG->vertices_num+2,sizeof(double));
    for(i=1;i<=TRG->vertices_num+1;i++) {
      TRG->cluster_degree[i]=SRC->cluster_degree[i];
    }
  }

  //if calc roles then init array or roles vector for all nodes
  if(GNRL_ST.calc_roles==TRUE && GNRL_ST.mtf_sz==3) {
    if( (TRG->roles_vec=(char**)calloc(TRG->vertices_num+1,sizeof(char*))) ==NULL)
      return RC_ERR;
    for(i=1;i<=TRG->vertices_num;i++)
      if( (TRG->roles_vec[i]=(char*)calloc(TOTAL_ROLES_3_NUM+1,sizeof(char))) ==NULL)
	return RC_ERR;
  }

  //copy hub params
  TRG->con_vertices_num=SRC->con_vertices_num;
  TRG->hub_deg=SRC->hub_deg;
  TRG->hub=SRC->hub;
  if(GNRL_ST.efc_prob){
    rc |= MatInit(&TRG->e_map,TRG->vertices_num,SPARSE);
    for(i=1; i<=TRG->edges_num; i++)
      MatAsgn(TRG->e_map,TRG->e_arr[i].s,TRG->e_arr[i].t,i);

    //if efc prob approach then allocate hub_edges_vec and hub_edges_indices
    TRG->hub_edges_vec=(int*)calloc(TRG->edges_num+1,sizeof(int));
    TRG->hub_edges_indices=(int*)calloc(TRG->hub_deg+1,sizeof(int));
    memcpy(TRG->hub_edges_vec,SRC->hub_edges_vec,(TRG->edges_num+1)*sizeof(int));
    memcpy(TRG->hub_edges_indices,SRC->hub_edges_indices,(TRG->hub_deg+1)*sizeof(int));

  }

  //if conserve layers mode then copy layers info
  if(GNRL_ST.r_conserve_layers==TRUE){
    TRG->num_of_layers=SRC->num_of_layers;
    TRG->layers_range=SRC->layers_range;

    TRG->layer_node_num=(int*)calloc(TRG->num_of_layers+1,sizeof(int));
    TRG->layer_edges_num=(int*)calloc(TRG->num_of_layers+1,sizeof(int));
    for(i=1;i<=TRG->num_of_layers;i++){
      TRG->layer_node_num[i]=SRC->layer_node_num[i];
      TRG->layer_edges_num[i]=SRC->layer_edges_num[i];
    }
  }

  //dump_network(TRG);
  *TRG_p=TRG;
  return RC_OK;
}

