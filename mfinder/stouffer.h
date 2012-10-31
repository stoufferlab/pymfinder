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
#ifndef __STOUFFER_H
#define __STOUFFER_H

#include "bits.h"
#include "common.h"

/******************* Definitions ****************************************/

// input for the mfinder wrapper
typedef struct {
  // in case we need to read in the network's links
  char* Filename;

  // in case we give the networks' link info as input
  int* Edges;
  unsigned int NumEdges;

  // for the motif_structure code
  // 1) what size motifs are we interested in?
  // 2) to how many randomizations should we compare the real network?
  unsigned int MotifSize;
  int NRandomizations;

  // for the participation code:
  // how many "motif members" can we keep track of?
  // to randomize the network beforehand
  // and to adjust the type of randomization
  int MaxMembersListSz;
  int Randomize;
  int UseMetropolis;
} mfinder_input;

/******************************* Functions *******************************/

void set_default_options();

int load_network_from_array(Network **N_p, int* edges, int edges_num);

int read_network(Network **N_p, mfinder_input mfinderi);

list64* motif_structure(mfinder_input mfinderi);

list64* motif_participation(mfinder_input mfinderi);

#endif
