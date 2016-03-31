/***************************************************************************
                       basics.h  -  global variables
                             -------------------
    begin                : 30-08-2004
    copyright            : (C) 2004 by Anke Busch
    email                : anke.busch@informatik.uni-freiburg.de
 ***************************************************************************/

#ifndef _BASICS__
#define _BASICS__

#include <ctype.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include "energy.h"
#include "interior_energy.h"
#include "stacking_energy.h"
#include "bulge_energy.h"
#include "dangling_energy.h"
#include "hybrid.h"
#include "single_run.h"

extern "C"{
             #include <ViennaRNA/fold.h>
             #include <ViennaRNA/part_func.h>
             #include <ViennaRNA/fold_vars.h>
             #include <ViennaRNA/utils.h>
             #include <ViennaRNA/LPfold.h>
             #include <ViennaRNA/PS_dot.h>
             #include <ViennaRNA/part_func_up.h>
          }

using namespace std;
/**********************************************************************************
*                    Structure- and Type definitions                              *
**********************************************************************************/
typedef unsigned int unt;


/**********************************************************************************
*                                 Makros                                          *
**********************************************************************************/
#define NUM_ACIDS 20
#define NUM_PAM 21


/**********************************************************************************
*                               Constants                                         *
**********************************************************************************/

extern const int MAX_BULGE_SIZE;
extern const int MAX_INTERIOR_ONE_SIDE_SIZE;

extern const double K_0;
extern const double T_MEASURE;
extern const double R;

extern const double MIN_DOUBLE;
extern const double MAX_DOUBLE;
extern const double DOUBLE_DIFF;

extern const int MIN_INT;
extern const int MAX_INT;


/**********************************************************************************
*                         global variabels (extern)                              *
**********************************************************************************/

extern double stacking_energies[];
extern double stacking_enthalpies[];
extern double interior_loop_1_1_energy[];
extern double interior_loop_1_1_enthalpy[];
extern double interior_loop_1_2_energy[];
extern double interior_loop_1_2_enthalpy[];
extern double interior_loop_2_2_energy[];
extern double interior_loop_2_2_enthalpy[];
extern double loop_destabilizing_energies[];
extern double single_base_stacking_energy[];
extern double single_base_stacking_enthalpies[];
extern double mismatch_energies_hairpin[];
extern double mismatch_energies_interior[];
extern double mismatch_enthalpies[];

extern double   RT                  ;
extern double   duplex_init         ;
extern double   terminal_AU         ;
extern double   loop_extrapolation  ;
extern double   max_ninio_correction;
extern double   ninio_correction    ;

extern char*    sequence_target     ;
extern char*    sequence_miRNA      ;
extern int*     int_sequence_target ;
extern int*     int_sequence_miRNA  ;
extern vector< string > targets     ;  
extern vector< string > miRNAs      ; 

extern int**    che1                ;
extern int**    che2                ;
extern int**    ches1               ;
extern int**    ches2               ;
extern double** C                   ;
extern double** Cseed               ;
extern double** U                   ;
extern double** ED_target           ;
extern bool     ED_target_computed  ;
extern double   ED_target_scaling   ;
extern double** ED_miRNA            ;
extern bool     ED_miRNA_computed   ;
extern double   ED_miRNA_scaling    ;
extern double** seed                ;
extern int**    seed_unpaired_target;
extern int**    seed_unpaired_miRNA ;
extern int      num_bp_seed         ;
extern int      max_unpaired_target ;
extern int      max_unpaired_miRNA  ;
extern int      max_unpaired        ;
// next line is quick hack, please remove if not needed anymore
extern double   pu_comb_threshold   ;
extern int      seed_region_start_ncRNA_5;
extern int      seed_region_end_ncRNA_5;
extern int      seed_region_start_ncRNA_3;
extern int      seed_region_end_ncRNA_3;
extern bool     seed_region_ncRNA   ;
extern bool     window              ;
extern int      window_size         ;
extern double   threshold           ;
extern bool     max_length          ;
extern int      max_unpaired_length ;
extern bool     max_dist            ;
extern int      max_pair_dist       ;
extern int      max_subopt          ;
extern int      best_che1           ;
extern int      best_che2           ;
extern bool     detailed_output     ;
extern bool     use_RNAplfold       ;
extern bool     use_RNAup           ;
extern bool     target_file         ;
extern bool     miRNA_file          ;
extern bool     heuristic           ;


extern double***** hybrid_matrix    ;

extern string   target_name         ; 
extern string   miRNA_name          ; 
extern string   hybrid_sequence_target;
extern string   hybrid_sequence_miRNA;
extern string   bulges_sequence_target;
extern string   bulges_sequence_miRNA;
extern string   hybrid_pos_target   ;

extern int      target_counter      ;
extern int      miRNA_counter       ;

extern int      pos_target_1        ;
extern int      pos_target_2        ;
extern int      pos_target_3        ;
extern int      pos_target_4        ;
extern int      pos_target_1_dangle ;
extern int      pos_target_4_dangle ;
extern int      pos_miRNA_1         ;
extern int      pos_miRNA_2         ;
extern int      pos_miRNA_3         ;
extern int      pos_miRNA_4         ;
extern int      pos_miRNA_1_dangle  ;
extern int      pos_miRNA_4_dangle  ;

/**********************************************************************************
*                         Funktionen                                              *
**********************************************************************************/

char int2char(int base);
char* int2char(int* num_seq, int size);
int char2int_base(char base);
int* char2int(char* seq);
void change_direction(char* sequence);

int BP2int(int b1, int b2);
void BP2_2(int bp, int &base_i, int &base_j);

bool PairTest(int b1, int b2);
bool PairTest(char b1, char b2);

bool check_sequence_and_set2upper(char* seq);
bool check_sequence_and_set2upper(string& seq);

bool is_ambiguous_sequence(char* seq);

double Double_Round(double x);
int Minimum(int a, int b);
double Minimum(double a, double b);
int Maximum(int a, int b);

double Sum_or_MAX_DOUBLE(double a, double b);
int SumVec(int* vec, int size);

double* MiniVec(double* vec, int size);
double* MaxiVec(double* vec, int size);


#endif   // _BASICS_
