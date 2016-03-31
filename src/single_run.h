#ifndef _SINGLE_RUN__
#define _SINGLE_RUN__

#include <stdlib.h>
#include "basics.h"

using namespace std;

bool readin_fasta(char* filename, vector<string>& vec);
bool test_options_hybridization();
bool test_options_window();
void init_target();
void init_miRNA();
void free_target();
void free_miRNA();
void clear_hydrization_seq();
void print_input();
void print_matrices();
void print_hybridization(double v);
void calculation();

#endif   // _SINGLE_RUN__
