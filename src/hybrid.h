#ifndef _HYBRID__
#define _HYBRID__

#include <stdlib.h>
#include "basics.h"

using namespace std;

void allocate_global_matrices();
void allocate_ED_target_matrix();
void allocate_ED_miRNA_matrix();
void init_global_matrices();
void init_ED_target_matrix();
void init_ED_miRNA_matrix();
void free_ED_target_matrix();
void free_ED_miRNA_matrix();
void free_global_matrices();

double calc_hybrid_value(int left_end_target, int left_end_miRNA, int num_bp, int unpaired_target, int unpaired_miRNA);
double compute_seed_ij(int i, int j);
double Cij     (int,int);
double Cij_seed(int,int);
double U_ij    (int,int);
double compute_hybrid() ;
pair<bool, double>  find_next_subopt();
void   traceback_right (int,int);
void   traceback_left  (int,int);
void   traceback_U     (int,int);
void   traceback_hybrid(int i, int j, int num_bp, int unpaired_target, int unpaired_miRNA);
void   traceback_seed  (int i, int j);
bool   traceback       (double);
void   traceback_print_labels_blanks(int i, int j);
void   traceback_print_hybrid_end(int i, int j);
#endif   // _HYBRID_
