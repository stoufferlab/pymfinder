/* mfinder.i */

/******************************* SWIGS **********************************/

%include <carrays.i>
%array_class(int, intArray);

%module mfinder
%{
#include "common.h"
#include "results.h"
#include "stouffer.h"
  extern list64* motif_structure(mfinder_input mfinderi);
%}

/******************************* Externs ********************************/

extern list64* motif_structure(mfinder_input mfinderi);

/******************************* Typedefs *******************************/

typedef struct {
  char* Filename;

  int* Edges;
  unsigned int NumEdges;

  unsigned int MotifSize;
  int NRandomizations;
} mfinder_input;

typedef long long int64;

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

/******************************* Functions *******************************/

%inline %{
    void motif_members(void* motif,
		       int* values,
		       int mtf_sz){
      unsigned i;
      Member* members = (Member*) motif;
      for(i=0;i<mtf_sz;++i){
	values[i] = members[i].node;
      }
    }

    Motif* get_motif(void* motif){
      return (Motif*) motif;
    }

    Motif_res* get_motif_res(void* motif_res){
      return (Motif_res*) motif_res;
    }
%}
