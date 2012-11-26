/* mfinder.i */

/******************************* SWIGS **********************************/

%include <carrays.i>
%array_class(int, intArray);

%module mfinder
%{
#include "common.h"
#include "results.h"
#include "stouffer.h"
  extern list* random_network(mfinder_input mfinderi);
  extern list64* motif_structure(mfinder_input mfinderi);
  extern list64* motif_participation(mfinder_input mfinderi);
%}

/******************************* Externs ********************************/

extern list* random_network(mfinder_input mfinderi);
extern list64* motif_structure(mfinder_input mfinderi);
extern list64* motif_participation(mfinder_input mfinderi);

/************************ Stouffer Typedefs *******************************/

typedef struct {
  char* Filename;
  int* Edges;
  unsigned int NumEdges;
  unsigned int MotifSize;
  int NRandomizations;
  int MaxMembersListSz;
  int Randomize;
  int UseMetropolis;
} mfinder_input;

/******************************* Typedefs *******************************/

typedef long long int64;

typedef struct _list_item{
	int val;
	void *p;
	struct _list_item *next;
} list_item;

typedef struct{
	int size;
	list_item *l;
} list;

typedef struct _list64_item{
	int64 val;
	void *p;
	struct _list64_item *next;
} list64_item;

typedef struct{
	int size;
	list64_item *l;
} list64;

// Statistics results about each subgraph
typedef struct {
	int64 id;
	int size;
	double real_count;
	double real_pval;
	double real_zscore;
	int hits_num;
	double conc_real;
	double conc_real_pval;
	double conc_real_zscore;
	//Neff_st *eff_arr;
	int64 unique_appear;
	double rand_mean;
	double rand_std_dev;
	double conc_rand_mean;
	double conc_rand_std_dev;
	list *all_members;
	int *all_rand_counts;
	double conv_grade;
	int dangling;
} Motif_res;

//Member info structure
typedef struct {
	int node;
	unsigned char role;
} Member;

//Motif info (Subgraph info) structure
typedef struct {
	int64 id;
	double count;
	double prob_count;
	double conc;
	int hits;
	list *members;
	list *all_members;
	int numberOfSelfEdges;
	double conv_grade;
} Motif;

// Edge structure
typedef struct {
	int s;
	int t;
	int weight;
} Edge;

/******************************* Functions *******************************/

%inline %{
    Edge* get_edge(void* p){
      return (Edge*) p;
    }

    Motif_res* get_motif_result(void* p){
      return (Motif_res*) p;
    }
    
    Motif* get_motif(void* p){
      return (Motif*) p;
    }

    void get_motif_members(void* p,
			   int* vals,
			   int mtf_sz){
      unsigned int i;
      Member* members = (Member*) p;
      for(i=0;i<mtf_sz;++i){
	vals[i] = members[i].node;
      }
    }

%}
