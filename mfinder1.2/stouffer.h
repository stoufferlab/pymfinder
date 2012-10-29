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

#include "bits.h"
#include "common.h"

/******************* Definitions ************************/

// input for the mfinder wrapper
typedef struct {
  char* Filename;

  int* Edges;
  unsigned int NumEdges;
  //int Weighted;

  unsigned int MotifSize;
  int NRandomizations;
} mfinder_input;

/******************************* Functions *******************************/

void set_default_options();
list64* motif_structure(mfinder_input mfinderi);
